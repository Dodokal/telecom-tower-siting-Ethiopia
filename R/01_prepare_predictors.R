# ============================================================================
# Ethiopia Telecom Tower Siting - Predictor Preparation (Non-GEE layers)
# ----------------------------------------------------------------------------
# Inputs (Windows paths, forward-slash-safe):
#   K:/ETH TOWERS/ETH_towers/                      <- GEE-derived rasters (template)
#   K:/ETH TOWERS/internet_users_norm_tifff/...    <- DHS mobile/internet users
#   K:/ETH TOWERS/electricity_Acess_norm/...       <- DHS electricity access
#   K:/ETH TOWERS/Ethiopia_flood_inundation/...    <- Trigg et al. flood inundation
#   K:/ETH TOWERS/ethiopia-260419-free/            <- Geofabrik OSM Ethiopia
#   K:/ETH TOWERS/grid.gpkg                        <- Arderne predictive MV grid (global)
#   K:/ETH TOWERS/towers_ethiopia_merged.gpkg      <- our merged tower presence set
#   K:/ETH TOWERS/2024-10-01_performance_mobile_tiles/gps_mobile_tiles.shp
#
# Outputs (written into K:/ETH TOWERS/ETH_towers/):
#   ETH_dhs_internet_users.tif
#   ETH_dhs_electricity_access.tif
#   ETH_flood_binary.tif
#   ETH_dist_to_road.tif
#   ETH_dist_to_primary_road.tif
#   ETH_dist_to_urban.tif
#   ETH_dist_to_airport.tif
#   ETH_dist_to_grid.tif
#   ETH_dist_to_existing_tower.tif
#   ETH_ookla_mobile_dl_kbps.tif
#
# Notes:
#  * All outputs share the GEE template grid (CRS, extent, resolution, alignment).
#  * Distance rasters are in metres (CRS is projected, EPSG:32637 expected).
#  * Continuous variables resampled with bilinear; categorical/binary with near.
#  * Designed to be re-runnable: each section can be commented in/out.
#
# Run-time on a normal laptop: roughly 20-40 minutes total. The slowest steps
# are the OSM road distance and Arderne grid clip.
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------------
required <- c("terra", "sf", "dplyr", "fs", "tibble")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(fs)
  library(tibble)
})

# Avoid GDAL/PROJ warnings cluttering the log
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 1, memfrac = 0.6)

# ---------------------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------------------
ROOT       <- "K:/ETH TOWERS"
GEE_DIR    <- file.path(ROOT, "ETH_towers")
OUT_DIR    <- GEE_DIR

DHS_INTERNET <- file.path(ROOT, "internet_users_norm_tifff",
                          "internet_users_norm_tifff_P6C7EhQ.tif")
DHS_ELEC     <- file.path(ROOT, "electricity_Acess_norm",
                          "electricity_Acess_norm_IrHaNqi.tif")
FLOOD_TIF    <- file.path(ROOT, "Ethiopia_flood_inundation",
                          "Ethiopia_flood_inundation_wcxTMEB.tif")

OSM_DIR        <- file.path(ROOT, "ethiopia-260419-free")
OSM_ROADS      <- file.path(OSM_DIR, "gis_osm_roads_free_1.shp")
OSM_PLACES     <- file.path(OSM_DIR, "gis_osm_places_free_1.shp")
OSM_TRANSPORT  <- file.path(OSM_DIR, "gis_osm_transport_a_free_1.shp")

GRID_GPKG    <- file.path(ROOT, "grid.gpkg")
TOWERS_GPKG  <- file.path(ROOT, "towers_ethiopia_merged.gpkg")
OOKLA_MOBILE <- file.path(ROOT, "2024-10-01_performance_mobile_tiles",
                          "gps_mobile_tiles.shp")

dir_create(OUT_DIR)

# ---------------------------------------------------------------------------
# 2. Load template raster (defines CRS / extent / resolution for everything)
# ---------------------------------------------------------------------------
pick_template <- function(dir) {
  preferred <- file.path(dir, "ETH_population.tif")
  if (file_exists(preferred)) return(preferred)
  cand <- dir_ls(dir, regexp = "\\.tif$", type = "file")
  if (length(cand) == 0)
    stop("No .tif files found in template directory: ", dir)
  as.character(cand[[1]])
}

TEMPLATE_PATH <- pick_template(GEE_DIR)
template <- rast(TEMPLATE_PATH)
template_crs <- crs(template)

cat("Template raster: ", TEMPLATE_PATH, "\n")
cat("  CRS    : ", crs(template, describe = TRUE)$name, "\n")
cat("  Extent : ", paste(round(as.vector(ext(template)), 0), collapse = ", "), "\n")
cat("  Res    : ", paste(round(res(template), 1), collapse = " x "), " (",
    crs(template, describe = TRUE)$name, ")\n", sep = "")
cat("  NCells : ", ncell(template), "\n")

# Helper: reproject + resample any raster onto the template grid
align_to_template <- function(r, method = "bilinear") {
  if (!same.crs(r, template))
    r <- project(r, template_crs, method = method, threads = TRUE)
  resample(r, template, method = method, threads = TRUE)
}

# Helper: build distance-to-feature raster (metres) from an sf vector
# Steps: reproject -> rasterize to template (1 where feature, NA elsewhere)
#        -> terra::distance returns metres to nearest non-NA cell
build_distance_raster <- function(vec, name, out_path) {
  cat("[", name, "] reproject + rasterize ...\n", sep = "")
  vec_proj <- st_transform(vec, st_crs(template_crs))
  vect_proj <- vect(vec_proj)
  r_burn <- rasterize(vect_proj, template, field = 1, background = NA,
                      touches = TRUE)
  cat("[", name, "] computing Euclidean distance ...\n", sep = "")
  d <- distance(r_burn)              # metres because CRS is projected
  d <- mask(d, template)             # keep extent identical to template
  names(d) <- name
  writeRaster(d, out_path, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))
  cat("[", name, "] -> ", out_path, "\n", sep = "")
  invisible(d)
}

# ---------------------------------------------------------------------------
# 3. DHS rasters (mobile/internet users + electricity access)
# ---------------------------------------------------------------------------
cat("\n=== 3. DHS rasters ===\n")

dhs_int  <- rast(DHS_INTERNET)
dhs_int  <- align_to_template(dhs_int, method = "bilinear")
names(dhs_int) <- "dhs_internet_users"
writeRaster(dhs_int,
            file.path(OUT_DIR, "ETH_dhs_internet_users.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

dhs_elec <- rast(DHS_ELEC)
dhs_elec <- align_to_template(dhs_elec, method = "bilinear")
names(dhs_elec) <- "dhs_electricity_access"
writeRaster(dhs_elec,
            file.path(OUT_DIR, "ETH_dhs_electricity_access.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# ---------------------------------------------------------------------------
# 4. Flood inundation -> binary (1 = flooded, 0 = not)
# ---------------------------------------------------------------------------
cat("\n=== 4. Flood inundation -> binary ===\n")
flood_raw <- rast(FLOOD_TIF)
flood_proj <- if (!same.crs(flood_raw, template))
  project(flood_raw, template_crs, method = "near", threads = TRUE) else flood_raw

# Anything > 0 (any inundation depth/extent) -> 1, else 0
flood_bin <- classify(flood_proj,
                      rcl = matrix(c(-Inf, 0,        0,
                                       0, Inf,       1), ncol = 3, byrow = TRUE),
                      include.lowest = TRUE, right = FALSE)
flood_bin <- subst(flood_bin, NA, 0)
flood_bin <- resample(flood_bin, template, method = "near")
flood_bin <- mask(flood_bin, template)
names(flood_bin) <- "flood_binary"
writeRaster(flood_bin,
            file.path(OUT_DIR, "ETH_flood_binary.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# ---------------------------------------------------------------------------
# 5. OSM roads -> dist_to_road (all) and dist_to_primary_road
# ---------------------------------------------------------------------------
cat("\n=== 5. OSM roads ===\n")
roads <- st_read(OSM_ROADS, quiet = TRUE)
cat("Roads loaded: ", nrow(roads), " features\n", sep = "")
print(table(roads$fclass)[1:min(15, length(unique(roads$fclass)))])

# All roads
build_distance_raster(
  vec      = roads,
  name     = "dist_to_road",
  out_path = file.path(OUT_DIR, "ETH_dist_to_road.tif"))

# Primary roads only: motorway / trunk / primary (and *_link variants)
primary_classes <- c("motorway", "motorway_link",
                     "trunk", "trunk_link",
                     "primary", "primary_link")
roads_primary <- roads %>% filter(fclass %in% primary_classes)
cat("Primary-class roads: ", nrow(roads_primary), " features\n", sep = "")

build_distance_raster(
  vec      = roads_primary,
  name     = "dist_to_primary_road",
  out_path = file.path(OUT_DIR, "ETH_dist_to_primary_road.tif"))

rm(roads, roads_primary); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 6. OSM places -> dist_to_urban (city + town centres)
# ---------------------------------------------------------------------------
cat("\n=== 6. OSM places (urban centres) ===\n")
places <- st_read(OSM_PLACES, quiet = TRUE)
print(sort(table(places$fclass), decreasing = TRUE))

urban <- places %>% filter(fclass %in% c("city", "town"))
cat("Urban centres (city+town): ", nrow(urban), "\n", sep = "")

build_distance_raster(
  vec      = urban,
  name     = "dist_to_urban",
  out_path = file.path(OUT_DIR, "ETH_dist_to_urban.tif"))

rm(places, urban); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 7. OSM transport (airports / aerodromes) -> dist_to_airport
# ---------------------------------------------------------------------------
cat("\n=== 7. OSM aerodromes ===\n")
transport <- st_read(OSM_TRANSPORT, quiet = TRUE)
print(sort(table(transport$fclass), decreasing = TRUE))

# Geofabrik tags: 'aerodrome' covers civilian airports; some extracts also
# expose 'airport' or 'helipad'. Take all aerodrome-like classes.
airport_classes <- intersect(unique(transport$fclass),
                             c("aerodrome", "airport", "airfield"))
if (length(airport_classes) == 0) airport_classes <- "aerodrome"

airports <- transport %>% filter(fclass %in% airport_classes)
cat("Aerodrome features: ", nrow(airports), "\n", sep = "")
if (nrow(airports) == 0)
  warning("No aerodrome features found. Skipping ETH_dist_to_airport.tif")

if (nrow(airports) > 0) {
  build_distance_raster(
    vec      = airports,
    name     = "dist_to_airport",
    out_path = file.path(OUT_DIR, "ETH_dist_to_airport.tif"))
}
rm(transport, airports); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 8. Arderne predictive MV grid (clip global -> Ethiopia, then distance)
# ---------------------------------------------------------------------------
cat("\n=== 8. Arderne predictive grid (clip + distance) ===\n")

# Build a generous Ethiopia bbox in WGS84 to clip the global grid efficiently
eth_bbox_4326 <- st_as_sfc(st_bbox(c(xmin = 32.5, ymin = 3.0,
                                     xmax = 48.5, ymax = 15.5),
                                   crs = 4326))
wkt_bbox <- st_as_text(eth_bbox_4326)

# Read only features intersecting the bbox (server-side filter via OGR)
grid_eth <- st_read(GRID_GPKG, wkt_filter = wkt_bbox, quiet = TRUE)
cat("Grid features in Ethiopia bbox: ", nrow(grid_eth), "\n", sep = "")

# Tighter clip just in case
grid_eth <- suppressWarnings(st_intersection(
  st_make_valid(grid_eth),
  st_transform(eth_bbox_4326, st_crs(grid_eth))))

build_distance_raster(
  vec      = grid_eth,
  name     = "dist_to_grid",
  out_path = file.path(OUT_DIR, "ETH_dist_to_grid.tif"))

rm(grid_eth); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 9. Existing towers -> dist_to_existing_tower (OSM-only subset)
# ---------------------------------------------------------------------------
cat("\n=== 9. Existing towers (OSM-only) ===\n")
towers_all <- st_read(TOWERS_GPKG, quiet = TRUE)
cat("Total tower records: ", nrow(towers_all), "\n", sep = "")
if ("sources" %in% names(towers_all)) {
  cat("By source:\n"); print(table(towers_all$sources))
  towers_osm <- towers_all %>%
    filter(grepl("OSM", sources, ignore.case = TRUE))
} else {
  message("No 'sources' column found - using all towers as-is.")
  towers_osm <- towers_all
}
cat("OSM-confirmed towers used: ", nrow(towers_osm), "\n", sep = "")

build_distance_raster(
  vec      = towers_osm,
  name     = "dist_to_existing_tower",
  out_path = file.path(OUT_DIR, "ETH_dist_to_existing_tower.tif"))

rm(towers_all, towers_osm); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 10. Ookla mobile performance (use mobile tiles, not fixed)
# ---------------------------------------------------------------------------
cat("\n=== 10. Ookla mobile performance ===\n")
# Reasoning: this study concerns mobile/cellular network siting.
# 'gps_fixed_tiles.shp' = wired broadband (DSL/fibre/cable) - irrelevant.
# 'gps_mobile_tiles.shp' = cellular (3G/4G/5G) - what we want.

ook_all <- st_read(OOKLA_MOBILE, quiet = TRUE)
cat("Ookla global mobile tiles: ", nrow(ook_all), "\n", sep = "")

# Clip to Ethiopia bbox
ook_eth <- suppressWarnings(
  st_intersection(st_make_valid(ook_all),
                  st_transform(eth_bbox_4326, st_crs(ook_all))))
cat("Ookla tiles within Ethiopia bbox: ", nrow(ook_eth), "\n", sep = "")
print(head(ook_eth))

# Field name can vary slightly between Ookla releases. Try the standard ones.
candidate_dl_fields <- c("avg_d_kbps", "avg_dl_kbps", "avg_dl",
                         "AvgDownloadKbps", "avg_download_kbps")
dl_field <- intersect(candidate_dl_fields, names(ook_eth))[1]
if (is.na(dl_field))
  stop("Could not find an average download field in the Ookla shapefile. ",
       "Inspect names(ook_eth) and update candidate_dl_fields.")
cat("Using download field: ", dl_field, "\n", sep = "")

ook_eth_proj <- st_transform(ook_eth, st_crs(template_crs))
ook_vect <- vect(ook_eth_proj)

ookla_r <- rasterize(ook_vect, template, field = dl_field, fun = "mean",
                     background = NA)
ookla_r <- mask(ookla_r, template)
names(ookla_r) <- "ookla_mobile_dl_kbps"
writeRaster(ookla_r,
            file.path(OUT_DIR, "ETH_ookla_mobile_dl_kbps.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

rm(ook_all, ook_eth, ook_eth_proj, ook_vect); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 11. Final inventory + sanity summary
# ---------------------------------------------------------------------------
cat("\n=== 11. Output inventory ===\n")
out_files <- dir_ls(OUT_DIR, regexp = "\\.tif$", type = "file")
inventory <- tibble(
  file       = path_file(out_files),
  size_MB    = round(as.numeric(file_info(out_files)$size) / 1024^2, 2),
  modified   = file_info(out_files)$modification_time
)
print(inventory)

write.csv(inventory,
          file.path(OUT_DIR, "_predictor_inventory.csv"),
          row.names = FALSE)

cat("\nAll predictor rasters written to: ", OUT_DIR, "\n", sep = "")
cat("Inventory CSV : ", file.path(OUT_DIR, "_predictor_inventory.csv"), "\n", sep = "")
cat("Done.\n")
