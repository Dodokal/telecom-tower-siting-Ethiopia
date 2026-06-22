# ============================================================================
# Ethiopia Telecom Tower Siting - RF + XGBoost suitability rasters
# ----------------------------------------------------------------------------
# Purpose:
#   The main pipeline (v4) only saved the LightGBM national raster to keep
#   runtime manageable. This standalone script:
#     1. Re-fits RF and XGBoost on the full training table (fast, ~2-5 min)
#     2. Predicts each at 1 km across Ethiopia (~10-20 min per model)
#     3. Writes results/07_suitability_rf_1km.tif and 07_suitability_xgb_1km.tif
#     4. Produces a 3-panel comparison figure: RF | XGBoost | LightGBM
#
# Inputs (all already on disk):
#   K:/ETH TOWERS/training_table.csv
#   K:/ETH TOWERS/ETH_towers/                 (predictor rasters)
#   K:/ETH TOWERS/results/07_suitability_1km.tif   (LightGBM, from v4)
#
# Outputs:
#   K:/ETH TOWERS/results/07_suitability_rf_1km.tif
#   K:/ETH TOWERS/results/07_suitability_xgb_1km.tif
#   K:/ETH TOWERS/results/figures/fig9_model_comparison.png
#
# Runtime: roughly 25-50 minutes total.
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------------
required <- c("terra", "sf", "dplyr", "tibble", "readr", "tidyr", "fs",
              "tidymodels", "ranger", "xgboost", "bonsai", "lightgbm",
              "ggplot2", "patchwork", "tidyterra", "viridis")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(fs)
  library(tidymodels); library(ranger); library(xgboost)
  library(bonsai); library(lightgbm)
  library(ggplot2); library(patchwork); library(tidyterra); library(viridis)
})

set.seed(42)
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 1, memfrac = 0.6)

# ---------------------------------------------------------------------------
# 1. Paths and config (must match the main pipeline)
# ---------------------------------------------------------------------------
ROOT       <- "K:/ETH TOWERS"
GEE_DIR    <- file.path(ROOT, "ETH_towers")
TRAIN_CSV  <- file.path(ROOT, "training_table.csv")
RES_DIR    <- file.path(ROOT, "results")
FIG_DIR    <- file.path(RES_DIR, "figures")
dir_create(RES_DIR); dir_create(FIG_DIR)

EXCLUDE_FROM_PREDICTORS <- c("ETH_ookla_mobile_dl_kbps",
                             "ETH_water_occurrence",
                             "ETH_dist_to_existing_tower")
EXCLUDE_FROM_RASTER_STACK <- c("ETH_ookla_mobile_dl_kbps",
                               "ETH_water_occurrence")
CATEGORICAL_NAMES <- c("ETH_landcover_2021", "ETH_smod_2020",
                       "ETH_protected_areas", "ETH_flood_binary")

PRED_AGGREGATION_FACT <- 10

theme_pub <- function() {
  theme_bw(base_size = 11, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"))
}
save_fig <- function(p, file, w, h) ggsave(file, p, width = w, height = h,
                                           dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# 2. Load training table + prep predictors
# ---------------------------------------------------------------------------
cat("=== 1. Load training table ===\n")
train_raw <- read_csv(TRAIN_CSV, show_col_types = FALSE)

META_COLS <- intersect(c("point_id", "label", "source",
                         "x_proj", "y_proj", "lon", "lat"),
                       names(train_raw))
predictor_cols <- setdiff(names(train_raw), META_COLS)
predictor_cols <- setdiff(predictor_cols, EXCLUDE_FROM_PREDICTORS)
train <- train_raw %>% drop_na(all_of(predictor_cols))

cat_cols  <- intersect(CATEGORICAL_NAMES, predictor_cols)
cont_cols <- setdiff(predictor_cols, cat_cols)
for (col in cat_cols) train[[col]] <- as.factor(train[[col]])

train$label_fct <- factor(ifelse(train$label == 1, "presence", "absence"),
                          levels = c("presence", "absence"))
final_predictors <- predictor_cols
cat("Training rows: ", nrow(train), ", predictors: ",
    length(final_predictors), "\n", sep = "")

# ---------------------------------------------------------------------------
# 3. Re-fit RF and XGBoost on full training set
# ---------------------------------------------------------------------------
cat("\n=== 2. Re-fit RF and XGBoost on full training set ===\n")
formula_str <- paste("label_fct ~", paste(final_predictors, collapse = " + "))
rec <- recipe(as.formula(formula_str), data = train) %>%
  step_novel(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

# Use sensible default hyperparameters (close to typical best from tuning).
# This avoids re-running the slow tune_grid step.
n_pred  <- length(final_predictors)
mtry_def <- max(2, floor(sqrt(n_pred)))

spec_rf <- rand_forest(trees = 1000, mtry = mtry_def, min_n = 5) %>%
  set_engine("ranger", num.threads = 4) %>%
  set_mode("classification")

spec_xgb <- boost_tree(trees = 800, tree_depth = 6, learn_rate = 0.05,
                       min_n = 5, sample_size = 0.85, mtry = mtry_def,
                       loss_reduction = 0) %>%
  set_engine("xgboost", nthread = 4) %>%
  set_mode("classification")

wf_rf  <- workflow() %>% add_recipe(rec) %>% add_model(spec_rf)
wf_xgb <- workflow() %>% add_recipe(rec) %>% add_model(spec_xgb)

cat("  fitting Random Forest ...\n")
fit_rf  <- fit(wf_rf,  data = train)
cat("  fitting XGBoost ...\n")
fit_xgb <- fit(wf_xgb, data = train)

# ---------------------------------------------------------------------------
# 4. Build 1km predictor stack
# ---------------------------------------------------------------------------
cat("\n=== 3. Build 1km predictor stack ===\n")
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
pred_stack_100m <- pred_stack_full[[intersect(final_predictors,
                                              names(pred_stack_full))]]

template_1km <- aggregate(pred_stack_100m[[1]],
                          fact = PRED_AGGREGATION_FACT,
                          fun = "mean", na.rm = TRUE)
cat("Resampling continuous predictors to 1km ...\n")
pred_stack_1km <- resample(pred_stack_100m, template_1km, method = "average",
                           threads = TRUE)

cat_layers_in_stack <- intersect(cat_cols, names(pred_stack_100m))
if (length(cat_layers_in_stack) > 0) {
  cat_resampled <- resample(pred_stack_100m[[cat_layers_in_stack]],
                            template_1km, method = "near", threads = TRUE)
  for (cl in cat_layers_in_stack)
    pred_stack_1km[[cl]] <- cat_resampled[[cl]]
}

country_mask <- pred_stack_1km[["ETH_population"]] >= 0
pred_stack_1km <- mask(pred_stack_1km, country_mask, maskvalue = 0)
cat("1km prediction stack ready: ", nlyr(pred_stack_1km), " layers, ",
    ncell(pred_stack_1km), " cells\n", sep = "")

# ---------------------------------------------------------------------------
# 5. Predict RF and XGBoost rasters
# ---------------------------------------------------------------------------
predict_chunk <- function(model, ...) {
  df <- as.data.frame(...)
  for (col in cat_cols) {
    if (col %in% names(df))
      df[[col]] <- factor(df[[col]], levels = levels(train[[col]]))
  }
  preds <- predict(model, new_data = df, type = "prob")
  preds$.pred_presence
}

cat("\n=== 4. Predict RF raster ===\n")
t0 <- Sys.time()
suit_rf_1km <- terra::predict(pred_stack_1km, fit_rf,
                              fun = predict_chunk, na.rm = TRUE,
                              cores = 1, progress = "text")
names(suit_rf_1km) <- "suitability_rf"
cat("RF prediction took ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " min\n", sep = "")
writeRaster(suit_rf_1km,
            file.path(RES_DIR, "07_suitability_rf_1km.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

cat("\n=== 5. Predict XGBoost raster ===\n")
t0 <- Sys.time()
suit_xgb_1km <- terra::predict(pred_stack_1km, fit_xgb,
                               fun = predict_chunk, na.rm = TRUE,
                               cores = 1, progress = "text")
names(suit_xgb_1km) <- "suitability_xgb"
cat("XGB prediction took ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    " min\n", sep = "")
writeRaster(suit_xgb_1km,
            file.path(RES_DIR, "07_suitability_xgb_1km.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# ---------------------------------------------------------------------------
# 6. Three-panel comparison figure (RF | XGBoost | LightGBM)
# ---------------------------------------------------------------------------
cat("\n=== 6. Build comparison figure ===\n")

# Load LightGBM raster from the main pipeline's output
lgbm_path <- file.path(RES_DIR, "07_suitability_1km.tif")
if (file_exists(lgbm_path)) {
  suit_lgbm_1km <- rast(lgbm_path)
  names(suit_lgbm_1km) <- "suitability_lgbm"
  has_lgbm <- TRUE
} else {
  message("LightGBM raster not found at ", lgbm_path,
          " - producing 2-panel figure instead.")
  has_lgbm <- FALSE
}

mk_panel <- function(r, title) {
  ggplot() +
    geom_spatraster(data = r) +
    scale_fill_viridis_c(name = "P(tower)", option = "magma",
                         na.value = "transparent", limits = c(0, 1),
                         guide = guide_colourbar(barwidth = unit(7, "cm"),
                                                 barheight = unit(0.3, "cm"))) +
    labs(title = title) +
    coord_sf() +
    theme_pub() +
    theme(panel.grid = element_blank(),
          axis.text  = element_text(size = 7),
          plot.title = element_text(face = "bold", size = 11))
}

p_rf  <- mk_panel(suit_rf_1km,  "Random Forest")
p_xgb <- mk_panel(suit_xgb_1km, "XGBoost")

if (has_lgbm) {
  p_lgbm <- mk_panel(suit_lgbm_1km, "LightGBM")
  combined <- (p_rf | p_xgb | p_lgbm) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Predicted telecom tower suitability, Ethiopia",
      subtitle = "Three machine-learning models at 1 km native resolution",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(colour = "grey30")))
  W <- 14; H <- 6
} else {
  combined <- (p_rf | p_xgb) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Predicted telecom tower suitability, Ethiopia",
      subtitle = "Two machine-learning models at 1 km native resolution",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(colour = "grey30")))
  W <- 11; H <- 6.5
}
combined <- combined & theme(legend.position = "bottom")

save_fig(combined, file.path(FIG_DIR, "fig9_model_comparison.png"),
         w = W, h = H)
cat("Saved -> ", file.path(FIG_DIR, "fig9_model_comparison.png"), "\n",
    sep = "")
cat("\nDONE.\n")
