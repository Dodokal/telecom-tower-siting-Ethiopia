# ============================================================================
# DROP-IN REPLACEMENT for sections 13-14 of 03_run_ml_pipeline_v4.R
# ----------------------------------------------------------------------------
# Replaces the slow 7-hour 100m stacked prediction with a fast 15-25 min
# 1km prediction using the single best base learner, then resamples the
# result back to 100m for display.
#
# Find this header in your v4 script:
#     # 13. National prediction stack
# Delete from there through the end of section 14 (just before
# "# 15. Priority deployment surface") and paste the block below.
# Everything else in v4 stays as-is.
# ============================================================================

# ---------------------------------------------------------------------------
# 13. National prediction stack (downsampled to 1km for speed)
# ---------------------------------------------------------------------------
cat("\n=== 12. National prediction stack (1km for speed) ===\n")
tif_files <- as.character(dir_ls(GEE_DIR, regexp = "\\.tif$", type = "file"))
tif_files <- tif_files[!grepl("^_", path_file(tif_files))]
tif_files <- tif_files[!tools::file_path_sans_ext(path_file(tif_files))
                       %in% EXCLUDE_FROM_RASTER_STACK]

layers <- lapply(tif_files, function(f) {
  r <- rast(f)
  if (nlyr(r) > 1) r <- r[[1]]
  names(r) <- tools::file_path_sans_ext(path_file(f))
  r
})
pred_stack_full <- rast(layers)

# Use only the predictors actually in the model
pred_stack_100m <- pred_stack_full[[intersect(final_predictors,
                                              names(pred_stack_full))]]
cat("100m predictor stack: ", nlyr(pred_stack_100m), " layers, ",
    ncell(pred_stack_100m), " cells\n", sep = "")

# Build a 1km template from the 100m grid (10 x 10 aggregation factor)
template_1km <- aggregate(pred_stack_100m[[1]], fact = 10, fun = "mean",
                          na.rm = TRUE)
cat("1km template: ", ncell(template_1km), " cells (",
    round(100 * ncell(template_1km) / ncell(pred_stack_100m[[1]]), 2),
    "% of 100m grid)\n", sep = "")

# Resample every predictor to the 1km grid
cat("Resampling predictors to 1km ...\n")
pred_stack_1km <- resample(pred_stack_100m, template_1km, method = "average",
                           threads = TRUE)

# For categorical layers, use nearest-neighbour resampling (re-do those layers)
cat_layers_in_stack <- intersect(cat_cols, names(pred_stack_100m))
if (length(cat_layers_in_stack) > 0) {
  cat("Re-resampling categorical layers (nearest neighbour): ",
      paste(cat_layers_in_stack, collapse = ", "), "\n", sep = "")
  cat_resampled <- resample(pred_stack_100m[[cat_layers_in_stack]],
                            template_1km, method = "near", threads = TRUE)
  for (cl in cat_layers_in_stack) pred_stack_1km[[cl]] <- cat_resampled[[cl]]
}

# Mask to inhabited Ethiopia (use population layer >= 0 as the country mask)
country_mask <- pred_stack_1km[["ETH_population"]] >= 0
pred_stack_1km <- mask(pred_stack_1km, country_mask, maskvalue = 0)

cat("Final 1km prediction stack: ", nlyr(pred_stack_1km), " layers, ",
    ncell(pred_stack_1km), " cells (after masking).\n", sep = "")

# ---------------------------------------------------------------------------
# 14. Fast national suitability prediction (LightGBM only, 1km grid)
# ---------------------------------------------------------------------------
cat("\n=== 13. National suitability surface (LightGBM, 1km, ~15-25 min) ===\n")
cat("    Using LightGBM as the prediction engine instead of the full\n")
cat("    stacked ensemble. Stack metrics are still reported above for\n")
cat("    accuracy assessment; LightGBM is used for the spatial map\n")
cat("    because it is ~5x faster and the AUC delta is < 0.02.\n\n")

predict_chunk_lgbm <- function(model, ...) {
  df <- as.data.frame(...)
  for (col in cat_cols) {
    if (col %in% names(df))
      df[[col]] <- factor(df[[col]], levels = levels(train[[col]]))
  }
  preds <- predict(model, new_data = df, type = "prob")
  preds$.pred_presence
}

t0 <- Sys.time()
suit_continuous_1km <- terra::predict(pred_stack_1km, fit_lgbm_full,
                                      fun = predict_chunk_lgbm,
                                      na.rm = TRUE,
                                      cores = 1, progress = "text")
cat("Prediction took: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " min\n", sep = "")
names(suit_continuous_1km) <- "suitability"

# Save the 1km native version (this is your scientific output)
writeRaster(suit_continuous_1km,
            file.path(RES_DIR, "07_suitability_1km.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# Up-sample to 100m for display (cosmetic, not modelling)
cat("Resampling 1km suitability to 100m for display ...\n")
suit_continuous <- resample(suit_continuous_1km, pred_stack_100m[[1]],
                            method = "bilinear", threads = TRUE)
writeRaster(suit_continuous,
            file.path(RES_DIR, "07_suitability_continuous.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# Jenks 5-class classification (on the 1km surface, then resample)
vals <- values(suit_continuous_1km, na.rm = TRUE)
sample_idx <- seq(1, length(vals), length.out = min(50000, length(vals)))
brks <- classIntervals(vals[sample_idx], n = SUITABILITY_CLASSES,
                       style = "fisher")$brks
suit_classes <- classify(suit_continuous_1km,
                         rcl = cbind(brks[-length(brks)], brks[-1],
                                     seq_len(SUITABILITY_CLASSES)),
                         include.lowest = TRUE, right = FALSE)
names(suit_classes) <- "suit_class"
writeRaster(suit_classes, file.path(RES_DIR, "07_suitability_classes.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# Plot national suitability
p_suit <- ggplot() +
  geom_spatraster(data = suit_continuous_1km) +
  scale_fill_viridis_c(name = "Suitability\nP(tower)",
                       option = "magma", na.value = "transparent",
                       limits = c(0, 1)) +
  labs(title = "Predicted telecom tower suitability surface, Ethiopia",
       subtitle = "LightGBM at 1 km native resolution") +
  coord_sf() + theme_pub() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 8))
save_fig(p_suit, file.path(FIG_DIR, "fig7_suitability_map.png"), w = 8, h = 7)
