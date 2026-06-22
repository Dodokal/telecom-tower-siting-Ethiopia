# ============================================================================
# Ethiopia Telecom Tower Siting - Regenerate pseudo-absences (clean v2)
# ----------------------------------------------------------------------------
# Why this script exists:
#   The original pseudo-absences were sampled from Ethiopia's bounding box,
#   which includes large swaths of Somalia, Eritrea, Sudan, Kenya, Djibouti
#   and the Gulf of Aden. About half (52%) of those absence points fell
#   outside the extent of the GEE predictor rasters and had to be dropped,
#   wrecking the presence/absence balance (424:205 instead of 424:424).
#
# What this script does:
#   1. Uses ETH_population.tif (or any GEE raster) as the strict country mask.
#      Pixels inside the raster extent and not NA = inside Ethiopia.
#   2. Resamples enough random candidate points within that mask.
#   3. Applies the same 1 km buffer-from-presence rule as before.
#   4. Writes a fresh presence_absence_v2.gpkg / .csv at 1:1 ratio.
#
# Then re-run training_table_v2.R to get a clean ~848-row training table.
# ============================================================================

required <- c("terra", "sf", "dplyr", "fs", "tibble", "readr")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
  library(fs); library(tibble); library(readr)
})

set.seed(42)
sf::sf_use_s2(FALSE)

# ---------------------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------------------
ROOT     <- "K:/ETH TOWERS"
GEE_DIR  <- file.path(ROOT, "ETH_towers")
TOWERS   <- file.path(ROOT, "towers_ethiopia_merged.gpkg")
OUT_GPKG <- file.path(ROOT, "presence_absence.gpkg")     # OVERWRITES original
OUT_CSV  <- file.path(ROOT, "presence_absence.csv")
BACKUP_TAG <- format(Sys.time(), "%Y%m%d_%H%M")

# Back up the old files first
if (file_exists(OUT_GPKG)) {
  file_move(OUT_GPKG,
            file.path(ROOT, paste0("presence_absence_old_", BACKUP_TAG, ".gpkg")))
}
if (file_exists(OUT_CSV)) {
  file_move(OUT_CSV,
            file.path(ROOT, paste0("presence_absence_old_", BACKUP_TAG, ".csv")))
}
cat("Backed up old presence_absence files with tag '", BACKUP_TAG, "'\n", sep = "")

# ---------------------------------------------------------------------------
# 2. Load Ethiopia mask from GEE raster
# ---------------------------------------------------------------------------
cat("\n=== 1. Build Ethiopia mask ===\n")

mask_candidates <- c("ETH_population.tif", "ETH_elevation.tif",
                     "ETH_landcover_2021.tif")
mask_path <- NULL
for (cand in mask_candidates) {
  p <- file.path(GEE_DIR, cand)
  if (file_exists(p)) { mask_path <- p; break }
}
if (is.null(mask_path)) {
  tifs <- as.character(dir_ls(GEE_DIR, regexp = "\\.tif$", type = "file"))
  if (length(tifs) == 0) stop("No GEE rasters in ", GEE_DIR)
  mask_path <- tifs[1]
}
cat("Using mask raster: ", path_file(mask_path), "\n", sep = "")

mask_r <- rast(mask_path)
if (nlyr(mask_r) > 1) mask_r <- mask_r[[1]]
mask_crs <- crs(mask_r)
cat("Mask CRS    : ", crs(mask_r, describe = TRUE)$name, "\n", sep = "")
cat("Mask extent : ", paste(round(as.vector(ext(mask_r)), 0), collapse = ", "),
    "\n", sep = "")

# Binary mask: 1 inside Ethiopia (non-NA), NA outside
inside_mask <- !is.na(mask_r)
inside_mask <- mask(inside_mask, mask_r)   # propagate NAs

n_inside_pixels <- sum(values(inside_mask), na.rm = TRUE)
cat("Pixels inside Ethiopia: ", n_inside_pixels, "\n", sep = "")

# ---------------------------------------------------------------------------
# 3. Load presence points (OSM towers)
# ---------------------------------------------------------------------------
cat("\n=== 2. Load presence points ===\n")
towers <- st_read(TOWERS, quiet = TRUE)

# Use OSM-confirmed towers as presences (this is what we agreed)
if ("sources" %in% names(towers)) {
  presence <- towers %>% filter(grepl("OSM", sources, ignore.case = TRUE))
} else {
  presence <- towers
}
presence <- st_transform(presence, mask_crs)
cat("OSM-confirmed presences: ", nrow(presence), "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Sample pseudo-absences strictly inside the country mask
# ---------------------------------------------------------------------------
cat("\n=== 3. Sample pseudo-absences inside mask ===\n")

n_target <- nrow(presence)            # 1:1 ratio
buffer_m <- 1000                      # 1 km buffer
oversample_factor <- 4                # we'll keep oversampling until we hit target

# Use spatSample with na.rm=TRUE to draw points only from non-NA pixels
candidates_needed <- n_target * oversample_factor
absence_pts <- vect()                 # empty SpatVector to grow

iter <- 0
while (nrow(absence_pts) < n_target && iter < 6) {
  iter <- iter + 1
  cat("  iter ", iter, ": sampling ", candidates_needed, " candidates...\n",
      sep = "")
  cand <- spatSample(inside_mask, size = candidates_needed,
                     method = "random", na.rm = TRUE,
                     as.points = TRUE, exhaustive = TRUE)
  if (is.null(cand) || nrow(cand) == 0) {
    cat("    spatSample returned 0 points; retrying with smaller request.\n")
    candidates_needed <- floor(candidates_needed / 2)
    next
  }
  # Drop candidates within 1 km of any presence
  pres_v <- vect(presence)
  d <- terra::distance(cand, pres_v, pairwise = FALSE)
  min_d <- apply(d, 1, min)
  cand_far <- cand[min_d >= buffer_m, ]
  cat("    ", nrow(cand_far), " of ", nrow(cand),
      " passed the ", buffer_m, " m buffer\n", sep = "")
  absence_pts <- rbind(absence_pts, cand_far)
  candidates_needed <- max((n_target - nrow(absence_pts)) * 4, 200)
}

if (nrow(absence_pts) < n_target) {
  warning("Could only generate ", nrow(absence_pts), " pseudo-absences ",
          "inside Ethiopia mask after ", iter, " iterations.")
}
absence_pts <- absence_pts[seq_len(min(n_target, nrow(absence_pts))), ]
cat("Pseudo-absences kept: ", nrow(absence_pts), "\n", sep = "")

# ---------------------------------------------------------------------------
# 5. Combine and write outputs
# ---------------------------------------------------------------------------
cat("\n=== 4. Assemble final presence/absence dataset ===\n")

# Convert absence points to sf
absence_sf <- st_as_sf(absence_pts) %>%
  mutate(label = 0L, source = "pseudo_absence")

# Drop any extra columns from presence that we don't need
presence_keep <- presence %>%
  st_geometry() %>%
  st_sf(geometry = .) %>%
  mutate(label = 1L, source = "OSM")

# Make sure both have only the columns we want
absence_clean <- absence_sf %>% select(label, source, geometry)
presence_clean <- presence_keep %>% select(label, source, geometry)

combined <- bind_rows(presence_clean, absence_clean)
combined <- st_transform(combined, 4326)
combined$point_id <- sprintf("PT_%05d", seq_len(nrow(combined)))
combined$lon <- st_coordinates(combined)[, 1]
combined$lat <- st_coordinates(combined)[, 2]
combined <- combined %>% select(point_id, label, source, lon, lat, geometry)

# Write
st_write(combined, OUT_GPKG, delete_dsn = TRUE, quiet = TRUE)
write_csv(st_drop_geometry(combined), OUT_CSV)

cat("\nFinal dataset:\n")
cat("  Total points : ", nrow(combined), "\n", sep = "")
cat("  Presences    : ", sum(combined$label == 1), "\n", sep = "")
cat("  Absences     : ", sum(combined$label == 0), "\n", sep = "")
cat("\nWritten to:\n  ", OUT_GPKG, "\n  ", OUT_CSV, "\n", sep = "")

# ---------------------------------------------------------------------------
# 6. Verify by re-extracting one predictor at the new points
# ---------------------------------------------------------------------------
cat("\n=== 5. Verify (no points should fall outside the mask) ===\n")
vec <- vect(combined)
vec_proj <- terra::project(vec, mask_crs)
test_extract <- terra::extract(mask_r, vec_proj, ID = FALSE)
n_na <- sum(is.na(test_extract[[1]]))
cat("Points outside mask: ", n_na, " (must be 0)\n", sep = "")
if (n_na > 0)
  warning("Some points fell outside the Ethiopia mask. Investigate.")

cat("\nDONE. Now re-run script 02 (build training table v2).\n")
cat("Expected: ~", nrow(combined), " rows, balanced ",
    sum(combined$label == 1), ":", sum(combined$label == 0), "\n", sep = "")
