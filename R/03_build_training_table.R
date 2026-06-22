# ============================================================================
# Ethiopia Telecom Tower Siting - Build training table (v3, robust)
# ----------------------------------------------------------------------------
# v3 changes vs v2:
#   - Extracts predictor rasters ONE LAYER AT A TIME instead of in bulk.
#     This sidesteps a terra bilinear-extraction edge case where bulk
#     extracts across many layers can return a different number of rows
#     than input points.
#   - Adds explicit nrow() checks at every step so any mismatch fails loud.
#   - Same predictor exclusions and imputation policy as v2.
# ============================================================================

required <- c("terra", "sf", "dplyr", "fs", "tibble", "readr")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
  library(fs); library(tibble); library(readr)
})

sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.6)

# ---------------------------------------------------------------------------
# Paths and config
# ---------------------------------------------------------------------------
ROOT     <- "K:/ETH TOWERS"
GEE_DIR  <- file.path(ROOT, "ETH_towers")
PA_GPKG  <- file.path(ROOT, "presence_absence.gpkg")
PA_CSV   <- file.path(ROOT, "presence_absence.csv")
OUT_CSV       <- file.path(ROOT, "training_table.csv")
OUT_TXT       <- file.path(ROOT, "training_table_summary.txt")
OUT_OOKLA_CSV <- file.path(ROOT, "ookla_validation.csv")

USE_CSV <- !file_exists(PA_GPKG)

CATEGORICAL <- c("ETH_landcover_2021",
                 "ETH_smod_2020",
                 "ETH_protected_areas",
                 "ETH_flood_binary")

EXCLUDE_FROM_PREDICTORS <- c("ETH_ookla_mobile_dl_kbps",
                             "ETH_water_occurrence")
VALIDATION_LAYER <- "ETH_ookla_mobile_dl_kbps"

MAX_NAS_PER_ROW <- 4

# ---------------------------------------------------------------------------
# 1. Load points
# ---------------------------------------------------------------------------
cat("=== 1. Load points ===\n")
if (USE_CSV) {
  pa <- read_csv(PA_CSV, show_col_types = FALSE)
  pa_sf <- st_as_sf(pa, coords = c("lon", "lat"), crs = 4326)
} else {
  pa_sf <- st_read(PA_GPKG, quiet = TRUE)
}
N_POINTS <- nrow(pa_sf)
cat("Points loaded: ", N_POINTS, "\n", sep = "")
print(table(pa_sf$label, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 2. Build predictor list (do NOT load as a single stack)
# ---------------------------------------------------------------------------
cat("\n=== 2. Predictor inventory ===\n")
tif_files <- as.character(dir_ls(GEE_DIR, regexp = "\\.tif$", type = "file"))
tif_files <- tif_files[!grepl("^_", path_file(tif_files))]

# Build a tibble: layer_name -> file_path
layer_index <- tibble(
  name = tools::file_path_sans_ext(path_file(tif_files)),
  path = tif_files
) %>%
  filter(!name %in% EXCLUDE_FROM_PREDICTORS)

cat("Predictor layers: ", nrow(layer_index), "\n", sep = "")

# Reference CRS = first layer's CRS
ref_r <- rast(layer_index$path[1])
if (nlyr(ref_r) > 1) ref_r <- ref_r[[1]]
ref_crs <- crs(ref_r)
cat("Reference CRS: ", crs(ref_r, describe = TRUE)$name, "\n", sep = "")

# ---------------------------------------------------------------------------
# 3. Reproject points ONCE
# ---------------------------------------------------------------------------
cat("\n=== 3. Reproject points to predictor CRS ===\n")
pts_proj <- st_transform(pa_sf, ref_crs)
pts_vect <- vect(pts_proj)
stopifnot(nrow(pts_vect) == N_POINTS)
cat("pts_vect rows: ", nrow(pts_vect), " (expected ", N_POINTS, ")\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Per-layer extraction with explicit row-count checks
# ---------------------------------------------------------------------------
cat("\n=== 4. Extract values, layer by layer ===\n")

extract_one_layer <- function(name, path, method) {
  r <- rast(path)
  if (nlyr(r) > 1) r <- r[[1]]
  names(r) <- name
  vals <- terra::extract(r, pts_vect, method = method, ID = FALSE)
  v <- vals[[1]]
  if (length(v) != N_POINTS) {
    # Defensive fallback: use simple method if bilinear misbehaves
    warning("Bilinear extract on ", name,
            " returned ", length(v), " rows; falling back to 'simple'.")
    vals <- terra::extract(r, pts_vect, method = "simple", ID = FALSE)
    v <- vals[[1]]
  }
  if (length(v) != N_POINTS) {
    stop("Layer ", name, " extract returned ", length(v),
         " rows; expected ", N_POINTS, ". Aborting.")
  }
  v
}

# Initialise the predictor data frame with one column per layer
pred_df <- tibble(.rows = N_POINTS)
for (i in seq_len(nrow(layer_index))) {
  ln <- layer_index$name[i]
  lp <- layer_index$path[i]
  meth <- if (ln %in% CATEGORICAL) "simple" else "bilinear"
  cat(sprintf("  [%2d/%2d] %-32s  method=%s\n",
              i, nrow(layer_index), ln, meth))
  pred_df[[ln]] <- extract_one_layer(ln, lp, meth)
}

stopifnot(nrow(pred_df) == N_POINTS)
cat("All layers extracted. pred_df: ", nrow(pred_df), " x ",
    ncol(pred_df), "\n", sep = "")

# ---------------------------------------------------------------------------
# 5. Ookla as separate validation layer
# ---------------------------------------------------------------------------
ookla_path <- file.path(GEE_DIR, paste0(VALIDATION_LAYER, ".tif"))
ookla_vals <- NULL
if (file_exists(ookla_path)) {
  cat("\n=== 5. Extract Ookla as validation layer ===\n")
  ookla_vals <- extract_one_layer(VALIDATION_LAYER, ookla_path, "bilinear")
  cat("Ookla coverage: ",
      sum(!is.na(ookla_vals)), " / ", N_POINTS,
      " (", round(100 * mean(!is.na(ookla_vals)), 1), "%)\n", sep = "")
}

# ---------------------------------------------------------------------------
# 6. Assemble training table
# ---------------------------------------------------------------------------
cat("\n=== 6. Assemble training table ===\n")
keep_cols <- intersect(c("point_id", "label", "source"), names(pts_proj))
pts_df <- st_drop_geometry(pts_proj)[, keep_cols, drop = FALSE]
pts_df$x_proj <- st_coordinates(pts_proj)[, 1]
pts_df$y_proj <- st_coordinates(pts_proj)[, 2]
pa_4326 <- st_transform(pts_proj, 4326)
pts_df$lon <- st_coordinates(pa_4326)[, 1]
pts_df$lat <- st_coordinates(pa_4326)[, 2]

stopifnot(nrow(pts_df) == nrow(pred_df))

train <- bind_cols(pts_df, pred_df)
train$label <- as.integer(train$label)

predictor_cols <- names(pred_df)
cat("Training table: ", nrow(train), " rows x ", ncol(train),
    " cols (", length(predictor_cols), " predictors)\n", sep = "")

# ---------------------------------------------------------------------------
# 7. Missingness diagnostics + impute + drop
# ---------------------------------------------------------------------------
cat("\n=== 7. Missingness (pre-imputation) ===\n")
miss_per_pred <- sapply(train[, predictor_cols, drop = FALSE],
                        function(x) sum(is.na(x)))
miss_per_pred_df <- tibble(
  predictor   = names(miss_per_pred),
  n_missing   = as.integer(miss_per_pred),
  pct_missing = round(100 * miss_per_pred / nrow(train), 2)
) %>% arrange(desc(n_missing))
print(miss_per_pred_df, n = Inf)

miss_per_point <- rowSums(is.na(train[, predictor_cols, drop = FALSE]))
cat("\nPoints with 0 NAs : ", sum(miss_per_point == 0), "\n", sep = "")
cat("Points 1-", MAX_NAS_PER_ROW, " NAs: ",
    sum(miss_per_point > 0 & miss_per_point <= MAX_NAS_PER_ROW), "\n", sep = "")
cat("Points > ", MAX_NAS_PER_ROW, " NAs: ",
    sum(miss_per_point > MAX_NAS_PER_ROW), " (will drop)\n", sep = "")

cat("\n=== 8. Impute + drop ===\n")
keep_rows <- miss_per_point <= MAX_NAS_PER_ROW
n_dropped <- sum(!keep_rows)
train_kept <- train[keep_rows, , drop = FALSE]
cat("Dropped ", n_dropped, " rows. Remaining: ",
    nrow(train_kept), "\n", sep = "")
print(table(train_kept$label))

mode_value <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0) return(NA)
  tt <- sort(table(ux), decreasing = TRUE)
  as(names(tt)[1], class(x))
}

cat_cols  <- intersect(CATEGORICAL, predictor_cols)
cont_cols <- setdiff(predictor_cols, cat_cols)

for (col in cont_cols) {
  if (any(is.na(train_kept[[col]]))) {
    med <- median(train_kept[[col]], na.rm = TRUE)
    train_kept[[col]][is.na(train_kept[[col]])] <- med
  }
}
for (col in cat_cols) {
  if (any(is.na(train_kept[[col]]))) {
    mv <- mode_value(train_kept[[col]])
    train_kept[[col]][is.na(train_kept[[col]])] <- mv
  }
}

stopifnot(sum(is.na(train_kept[, predictor_cols])) == 0)

# ---------------------------------------------------------------------------
# 8. Write outputs
# ---------------------------------------------------------------------------
cat("\n=== 9. Write outputs ===\n")
write_csv(train_kept, OUT_CSV)
cat("Training table -> ", OUT_CSV, "\n", sep = "")

if (!is.null(ookla_vals)) {
  ookla_df <- pts_df %>% mutate(ookla_dl_kbps = ookla_vals) %>%
    filter(!is.na(ookla_dl_kbps))
  write_csv(ookla_df, OUT_OOKLA_CSV)
  cat("Ookla validation -> ", OUT_OOKLA_CSV,
      " (", nrow(ookla_df), " rows)\n", sep = "")
}

sink(OUT_TXT)
cat("Ethiopia Telecom Tower Siting - Training table summary (v3)\n")
cat("Generated: ", format(Sys.time()), "\n", sep = "")
cat("------------------------------------------------------\n\n")
cat("Final training table\n--------------------\n")
cat("Rows                    : ", nrow(train_kept), "\n", sep = "")
cat("Columns (total)         : ", ncol(train_kept), "\n", sep = "")
cat("Predictors (continuous) : ", length(cont_cols), "\n", sep = "")
cat("Predictors (categorical): ", length(cat_cols), "\n", sep = "")
cat("\nLabel distribution:\n"); print(table(train_kept$label))
cat("\nMissingness BEFORE imputation:\n")
print(miss_per_pred_df, n = Inf)
cat("\nDropped rows: ", n_dropped, "\n", sep = "")
sink()

cat("\nDONE. Final: ", nrow(train_kept), " rows, ",
    ncol(train_kept), " cols.\n", sep = "")
