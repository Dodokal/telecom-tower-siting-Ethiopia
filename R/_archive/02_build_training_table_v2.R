# ============================================================================
# Ethiopia Telecom Tower Siting - Build training table (REVISED v2)
# ----------------------------------------------------------------------------
# Changes vs v1:
#   - Drops ETH_ookla_mobile_dl_kbps from predictors (95.5% NA in v1 run)
#     and saves it separately as an EXTERNAL VALIDATION layer.
#   - Drops ETH_water_occurrence from predictors (sampled all zeros in v1).
#   - Median-imputes 1-3 NA values per row instead of deleting whole rows
#     so we keep the full presence/absence balance.
#
# Output:
#   K:/ETH TOWERS/training_table.csv          <- model training input
#   K:/ETH TOWERS/training_table_summary.txt  <- diagnostics
#   K:/ETH TOWERS/ookla_validation.csv        <- separate validation layer
# ============================================================================

required <- c("terra", "sf", "dplyr", "fs", "tibble", "readr")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
  library(fs); library(tibble); library(readr)
})

sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 1, memfrac = 0.6)

# ---------------------------------------------------------------------------
# 1. Paths and config
# ---------------------------------------------------------------------------
ROOT     <- "K:/ETH TOWERS"
GEE_DIR  <- file.path(ROOT, "ETH_towers")
PA_GPKG  <- file.path(ROOT, "presence_absence.gpkg")
PA_CSV   <- file.path(ROOT, "presence_absence.csv")
OUT_CSV       <- file.path(ROOT, "training_table.csv")
OUT_TXT       <- file.path(ROOT, "training_table_summary.txt")
OUT_OOKLA_CSV <- file.path(ROOT, "ookla_validation.csv")

USE_CSV <- !file_exists(PA_GPKG)

# Layers handled as categorical at extraction (use nearest-neighbour)
CATEGORICAL <- c("ETH_landcover_2021",
                 "ETH_smod_2020",
                 "ETH_protected_areas",
                 "ETH_flood_binary")

# Layers EXCLUDED from predictor table:
#   - Ookla: 95% missing -> external validation only
#   - Water occurrence: all zeros over Ethiopia in this dataset
EXCLUDE_FROM_PREDICTORS <- c("ETH_ookla_mobile_dl_kbps",
                             "ETH_water_occurrence")
# Layer kept separately as external validation layer
VALIDATION_LAYER <- "ETH_ookla_mobile_dl_kbps"

# Imputation policy
MAX_NAS_PER_ROW <- 4   # rows with > this many NAs are dropped (likely outside study area)

# ---------------------------------------------------------------------------
# 2. Load points
# ---------------------------------------------------------------------------
cat("=== 1. Load points ===\n")
if (USE_CSV) {
  pa <- read_csv(PA_CSV, show_col_types = FALSE)
  pa_sf <- st_as_sf(pa, coords = c("lon", "lat"), crs = 4326)
} else {
  pa_sf <- st_read(PA_GPKG, quiet = TRUE)
}
cat("Loaded points: ", nrow(pa_sf), "\n", sep = "")
cat("Label distribution:\n"); print(table(pa_sf$label, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 3. Build predictor stack (with multi-band fix)
# ---------------------------------------------------------------------------
cat("\n=== 2. Build predictor stack ===\n")
tif_files <- as.character(dir_ls(GEE_DIR, regexp = "\\.tif$", type = "file"))
tif_files <- tif_files[!grepl("^_", path_file(tif_files))]
if (length(tif_files) == 0) stop("No .tif files in ", GEE_DIR)

cat("Files found in folder: ", length(tif_files), "\n", sep = "")

# Load each file -> keep only first band -> name from filename
all_layers <- lapply(tif_files, function(f) {
  r <- rast(f)
  if (nlyr(r) > 1) r <- r[[1]]
  names(r) <- tools::file_path_sans_ext(path_file(f))
  r
})
all_stk <- rast(all_layers)

# Separate: predictor stack vs Ookla (validation only)
predictor_names <- setdiff(names(all_stk), EXCLUDE_FROM_PREDICTORS)
stk_pred <- all_stk[[predictor_names]]

cat("Predictor layers (used in model): ", nlyr(stk_pred), "\n", sep = "")
for (n in names(stk_pred)) cat("  ", n, "\n", sep = "")
cat("\nExcluded from predictors: ",
    paste(intersect(names(all_stk), EXCLUDE_FROM_PREDICTORS), collapse = ", "),
    "\n", sep = "")
cat("  (Ookla saved separately as external validation layer.)\n")

# ---------------------------------------------------------------------------
# 4. Reproject points
# ---------------------------------------------------------------------------
cat("\n=== 3. Reproject points to stack CRS ===\n")
pts_proj <- st_transform(pa_sf, crs(stk_pred))
pts_vect <- vect(pts_proj)

# ---------------------------------------------------------------------------
# 5. Extract predictor values (continuous + categorical)
# ---------------------------------------------------------------------------
cat("\n=== 4. Extract values ===\n")
cat_layers  <- intersect(CATEGORICAL, names(stk_pred))
cont_layers <- setdiff(names(stk_pred), cat_layers)
cat("Continuous : ", length(cont_layers), "\n", sep = "")
cat("Categorical: ", length(cat_layers),  "\n", sep = "")

ext_cont <- terra::extract(stk_pred[[cont_layers]], pts_vect,
                           method = "bilinear", ID = FALSE)
ext_cat  <- if (length(cat_layers) > 0)
  terra::extract(stk_pred[[cat_layers]], pts_vect,
                 method = "simple", ID = FALSE)
  else NULL

# ---------------------------------------------------------------------------
# 6. Extract Ookla SEPARATELY (validation layer, not a predictor)
# ---------------------------------------------------------------------------
if (VALIDATION_LAYER %in% names(all_stk)) {
  cat("\n=== 5. Extract Ookla as validation layer ===\n")
  ext_ookla <- terra::extract(all_stk[[VALIDATION_LAYER]], pts_vect,
                              method = "bilinear", ID = FALSE)
  cat("Ookla coverage at points: ",
      sum(!is.na(ext_ookla[[1]])), " / ", nrow(ext_ookla),
      " (", round(100 * mean(!is.na(ext_ookla[[1]])), 1),
      "%)\n", sep = "")
}

# ---------------------------------------------------------------------------
# 7. Assemble training table
# ---------------------------------------------------------------------------
cat("\n=== 6. Assemble training table ===\n")
keep_cols <- intersect(c("point_id", "label", "source"), names(pts_proj))
pts_df <- st_drop_geometry(pts_proj)[, keep_cols, drop = FALSE]
pts_df$x_proj <- st_coordinates(pts_proj)[, 1]
pts_df$y_proj <- st_coordinates(pts_proj)[, 2]
pa_4326 <- st_transform(pts_proj, 4326)
pts_df$lon <- st_coordinates(pa_4326)[, 1]
pts_df$lat <- st_coordinates(pa_4326)[, 2]

train <- bind_cols(pts_df, ext_cont)
if (!is.null(ext_cat)) train <- bind_cols(train, ext_cat)

train$label <- as.integer(train$label)

predictor_cols <- c(cont_layers, cat_layers)
cat("Training table: ", nrow(train), " rows x ", ncol(train),
    " cols (", length(predictor_cols), " predictors)\n", sep = "")

# ---------------------------------------------------------------------------
# 8. Diagnostics: missingness BEFORE imputation
# ---------------------------------------------------------------------------
cat("\n=== 7. Missingness diagnostics (pre-imputation) ===\n")
miss_per_pred <- sapply(train[, predictor_cols, drop = FALSE],
                        function(x) sum(is.na(x)))
miss_per_pred_df <- tibble(
  predictor   = names(miss_per_pred),
  n_missing   = as.integer(miss_per_pred),
  pct_missing = round(100 * miss_per_pred / nrow(train), 2)
) %>% arrange(desc(n_missing))
print(miss_per_pred_df, n = Inf)

miss_per_point <- rowSums(is.na(train[, predictor_cols, drop = FALSE]))
cat("\nPoints with 0 NAs       : ", sum(miss_per_point == 0), "\n", sep = "")
cat("Points with 1-",MAX_NAS_PER_ROW," NAs    : ",
    sum(miss_per_point > 0 & miss_per_point <= MAX_NAS_PER_ROW), "\n", sep = "")
cat("Points with > ", MAX_NAS_PER_ROW, " NAs   : ",
    sum(miss_per_point > MAX_NAS_PER_ROW), "  (will be DROPPED)\n", sep = "")

# ---------------------------------------------------------------------------
# 9. Drop rows with too many NAs, then median-impute the rest
# ---------------------------------------------------------------------------
cat("\n=== 8. Impute + drop ===\n")
keep_rows <- miss_per_point <= MAX_NAS_PER_ROW
n_dropped <- sum(!keep_rows)
train_kept <- train[keep_rows, , drop = FALSE]
cat("Dropped ", n_dropped, " rows.\n", sep = "")
cat("Remaining rows: ", nrow(train_kept), "\n", sep = "")
cat("Label balance after drop:\n"); print(table(train_kept$label))

# Median-impute remaining NAs (continuous) and mode-impute (categorical)
mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0) return(NA)
  tt <- sort(table(ux), decreasing = TRUE)
  as(names(tt)[1], class(x))
}

cat("\nImputing remaining NAs (median for continuous, mode for categorical):\n")
for (col in cont_layers) {
  if (any(is.na(train_kept[[col]]))) {
    med <- median(train_kept[[col]], na.rm = TRUE)
    n_imp <- sum(is.na(train_kept[[col]]))
    train_kept[[col]][is.na(train_kept[[col]])] <- med
    cat("  ", col, ": ", n_imp, " imputed (median = ",
        round(med, 2), ")\n", sep = "")
  }
}
for (col in cat_layers) {
  if (any(is.na(train_kept[[col]]))) {
    mv <- mode_value(train_kept[[col]])
    n_imp <- sum(is.na(train_kept[[col]]))
    train_kept[[col]][is.na(train_kept[[col]])] <- mv
    cat("  ", col, ": ", n_imp, " imputed (mode = ", mv, ")\n", sep = "")
  }
}

# Sanity: zero NAs left in predictors
remaining_na <- sum(is.na(train_kept[, predictor_cols]))
cat("\nNAs remaining in predictor columns: ", remaining_na,
    " (must be 0)\n", sep = "")
stopifnot(remaining_na == 0)

# ---------------------------------------------------------------------------
# 10. Write outputs
# ---------------------------------------------------------------------------
cat("\n=== 9. Write outputs ===\n")
write_csv(train_kept, OUT_CSV)
cat("Training table -> ", OUT_CSV, "\n", sep = "")

# Ookla validation table (kept separately, point_id-keyed)
if (exists("ext_ookla")) {
  ookla_df <- pts_df %>%
    mutate(ookla_dl_kbps = ext_ookla[[1]]) %>%
    filter(!is.na(ookla_dl_kbps))
  write_csv(ookla_df, OUT_OOKLA_CSV)
  cat("Ookla validation -> ", OUT_OOKLA_CSV,
      " (", nrow(ookla_df), " rows with Ookla data)\n", sep = "")
}

# Plain-text summary
sink(OUT_TXT)
cat("Ethiopia Telecom Tower Siting - Training table summary (v2)\n")
cat("Generated: ", format(Sys.time()), "\n", sep = "")
cat("------------------------------------------------------\n\n")
cat("Pipeline changes vs v1:\n")
cat("  - Dropped Ookla from predictors (95% NA in v1); saved as validation layer.\n")
cat("  - Dropped ETH_water_occurrence (sampled all zeros).\n")
cat("  - Median/mode-imputed remaining NAs instead of deleting rows.\n\n")
cat("Final training table\n")
cat("--------------------\n")
cat("Rows                  : ", nrow(train_kept), "\n", sep = "")
cat("Columns (total)       : ", ncol(train_kept), "\n", sep = "")
cat("Predictors (continuous): ", length(cont_layers), "\n", sep = "")
cat("Predictors (categorical): ", length(cat_layers), "\n", sep = "")
cat("\nLabel distribution:\n"); print(table(train_kept$label))
cat("\nMissingness BEFORE imputation:\n")
print(miss_per_pred_df, n = Inf)
cat("\nDropped rows (>", MAX_NAS_PER_ROW, " NAs): ", n_dropped, "\n", sep = "")
sink()
cat("Summary       -> ", OUT_TXT, "\n", sep = "")

cat("\nDone. Final table: ", nrow(train_kept), " rows, ",
    ncol(train_kept), " cols.\n", sep = "")
