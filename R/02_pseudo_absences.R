# ============================================================================
# Ethiopia Telecom Tower Siting - Pseudo-absences v4
# ----------------------------------------------------------------------------
# Why this script exists:
#   The previous pseudo-absences were sampled uniformly across Ethiopia.
#   Towers cluster in populated areas; uniform absences are mostly empty
#   bushland. This makes the classifier's job artificially easy - it just
#   learns "is this a city or not". Spatial models trained that way over-
#   estimate AUC dramatically and reviewers at Q1 journals will catch it.
#
# What changes in v4:
#   Sample pseudo-absences with probability proportional to log-transformed
#   population density. This is the "target-group background" approach from
#   Phillips et al. 2009 (Ecological Applications, 19(1), 181-197). It forces
#   the classifier to distinguish between LIKELY tower sites in populated
#   areas vs. ACTUAL tower sites - a much harder and more honest task.
#
# Same constraints as before:
#   - Strictly inside Ethiopia (raster mask)
#   - >=1000 m from any presence (avoids existing service footprints)
#   - 1:1 ratio with presences
#
# Output: K:/ETH TOWERS/presence_absence.gpkg + .csv  (overwrites previous)
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
# Paths
# ---------------------------------------------------------------------------
ROOT     <- "K:/ETH TOWERS"
GEE_DIR  <- file.path(ROOT, "ETH_towers")
TOWERS   <- file.path(ROOT, "towers_ethiopia_merged.gpkg")
POP_TIF  <- file.path(GEE_DIR, "ETH_population.tif")
OUT_GPKG <- file.path(ROOT, "presence_absence.gpkg")
OUT_CSV  <- file.path(ROOT, "presence_absence.csv")
BACKUP_TAG <- format(Sys.time(), "v3_%Y%m%d_%H%M")

# Back up previous files
if (file_exists(OUT_GPKG))
  file_move(OUT_GPKG,
            file.path(ROOT, paste0("presence_absence_", BACKUP_TAG, ".gpkg")))
if (file_exists(OUT_CSV))
  file_move(OUT_CSV,
            file.path(ROOT, paste0("presence_absence_", BACKUP_TAG, ".csv")))
cat("Backed up old files with tag '", BACKUP_TAG, "'\n", sep = "")

# ---------------------------------------------------------------------------
# 1. Load population raster as the sampling-weight surface
# ---------------------------------------------------------------------------
cat("\n=== 1. Build population-weighted sampling surface ===\n")
pop <- rast(POP_TIF)
if (nlyr(pop) > 1) pop <- pop[[1]]
mask_crs <- crs(pop)

# Convert to sampling weights: log1p(pop) so dense cities don't dominate
weights_r <- log1p(clamp(pop, 0, NA))

# Floor: every inhabited pixel gets at least a tiny baseline probability so
# the model also sees some genuinely rural absences (~5% of the weight mass).
nz_vals <- values(weights_r, na.rm = TRUE)
nz_vals <- nz_vals[nz_vals > 0]
floor_val <- quantile(nz_vals, 0.10)        # 10th percentile of inhabited cells
weights_r <- ifel(is.na(weights_r), NA,
                  pmax(weights_r, floor_val))
# Inside-country mask (population raster is NA outside Ethiopia)
weights_r <- mask(weights_r, pop)

cat("Sampling weights summary:\n")
print(summary(values(weights_r, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 2. Load presences (OSM-confirmed)
# ---------------------------------------------------------------------------
cat("\n=== 2. Load presences ===\n")
towers <- st_read(TOWERS, quiet = TRUE)
if ("sources" %in% names(towers)) {
  presence <- towers %>% filter(grepl("OSM", sources, ignore.case = TRUE))
} else {
  presence <- towers
}
presence <- st_transform(presence, mask_crs)
cat("OSM-confirmed presences: ", nrow(presence), "\n", sep = "")

# ---------------------------------------------------------------------------
# 3. Population-weighted sampling of pseudo-absences
# ---------------------------------------------------------------------------
cat("\n=== 3. Sample pseudo-absences (population-weighted) ===\n")

n_target <- nrow(presence)
buffer_m <- 1000
oversample <- 6

absence_pts <- vect()
iter <- 0
while (nrow(absence_pts) < n_target && iter < 6) {
  iter <- iter + 1
  to_draw <- max((n_target - nrow(absence_pts)) * oversample, 500)
  cat("  iter ", iter, ": drawing ", to_draw, " weighted candidates...\n",
      sep = "")
  cand <- spatSample(weights_r, size = to_draw,
                     method = "weights", na.rm = TRUE,
                     as.points = TRUE, exhaustive = TRUE)
  if (is.null(cand) || nrow(cand) == 0) {
    cat("    nothing returned; halving request and retrying.\n")
    oversample <- max(2, oversample - 1); next
  }
  # Buffer rule: drop candidates within 1000m of any presence
  pres_v <- vect(presence)
  d <- terra::distance(cand, pres_v, pairwise = FALSE)
  min_d <- apply(d, 1, min)
  cand_far <- cand[min_d >= buffer_m, ]
  cat("    ", nrow(cand_far), " of ", nrow(cand), " kept after buffer\n",
      sep = "")
  absence_pts <- rbind(absence_pts, cand_far)
}

absence_pts <- absence_pts[seq_len(min(n_target, nrow(absence_pts))), ]
cat("Pseudo-absences: ", nrow(absence_pts), "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Combine and write
# ---------------------------------------------------------------------------
cat("\n=== 4. Write presence_absence files ===\n")
absence_sf <- st_as_sf(absence_pts) %>%
  mutate(label = 0L, source = "pseudo_absence") %>%
  select(label, source, geometry)

presence_sf <- st_geometry(presence) %>% st_sf(geometry = .) %>%
  mutate(label = 1L, source = "OSM") %>%
  select(label, source, geometry)

combined <- bind_rows(presence_sf, absence_sf) %>%
  st_transform(4326)
combined$point_id <- sprintf("PT_%05d", seq_len(nrow(combined)))
combined$lon <- st_coordinates(combined)[, 1]
combined$lat <- st_coordinates(combined)[, 2]
combined <- combined %>% select(point_id, label, source, lon, lat, geometry)

st_write(combined, OUT_GPKG, delete_dsn = TRUE, quiet = TRUE)
write_csv(st_drop_geometry(combined), OUT_CSV)

cat("Final dataset:\n")
cat("  Total points: ", nrow(combined), "\n", sep = "")
cat("  Presences   : ", sum(combined$label == 1), "\n", sep = "")
cat("  Absences    : ", sum(combined$label == 0), "\n", sep = "")
cat("\nWritten to:\n  ", OUT_GPKG, "\n  ", OUT_CSV, "\n", sep = "")
cat("\nNow re-run the training-table builder (script 02 v3) on the new file.\n")
