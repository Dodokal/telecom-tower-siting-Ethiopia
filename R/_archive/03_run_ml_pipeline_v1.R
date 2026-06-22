# ============================================================================
# Ethiopia Telecom Tower Siting - Full ML pipeline
# ----------------------------------------------------------------------------
# Reads:
#   K:/ETH TOWERS/training_table.csv          (from script 02)
#   K:/ETH TOWERS/ETH_towers/                 (predictor rasters)
#
# Writes (all into K:/ETH TOWERS/results/):
#   01_vif_screening.csv
#   02_cv_metrics_per_fold.csv
#   03_cv_metrics_summary.csv
#   04_wilcoxon_pairwise.csv
#   05_var_importance_rf.csv
#   05_var_importance_xgb.csv
#   05_var_importance_lgbm.csv
#   06_shap_values.csv
#   07_suitability_continuous.tif       (national 0-1 surface)
#   07_suitability_classes.tif          (5-class Jenks)
#   08_priority_greenfield.tif          (binary)
#   08_priority_densification.tif       (binary)
#   08_priority_score.tif               (continuous priority index)
#   09_run_summary.txt
#
# What it runs:
#   3.  VIF + |rho|>=0.85 multicollinearity screening
#   4.  blockCV spatial folds (50 km, 10 outer)
#   5.  Random Forest, XGBoost, LightGBM, MaxEnt (tuned via spatial CV)
#   6.  Stacked ensemble (RF + XGB + LGBM via elastic-net meta-learner)
#   7.  Performance metrics: AUC, AUC-PR, TSS, Kappa, F1, sens, spec, bal-acc
#   8.  Wilcoxon paired tests vs the stacked ensemble
#   9.  Permutation + SHAP interpretation
#   10. Full-country prediction -> 100m suitability raster + 5 Jenks classes
#   11. Priority deployment surface: greenfield vs densification
#
# Runtime: 1-3 hours depending on machine. Tunable via TUNE_GRID_SIZE below.
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Packages (auto-install missing ones)
# ---------------------------------------------------------------------------
required <- c(
  "terra", "sf", "dplyr", "tibble", "readr", "tidyr", "purrr", "fs",
  "tidymodels", "ranger", "xgboost", "lightgbm", "bonsai",
  "stacks", "yardstick",
  "blockCV", "PresenceAbsence",
  "fastshap", "vip",
  "classInt"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(purrr); library(fs)
  library(tidymodels); library(ranger); library(xgboost)
  library(lightgbm); library(bonsai)            # bonsai = lightgbm engine for parsnip
  library(stacks); library(yardstick)
  library(blockCV); library(PresenceAbsence)
  library(fastshap); library(vip)
  library(classInt)
})

set.seed(42)
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 1, memfrac = 0.6)

# ---------------------------------------------------------------------------
# 1. Paths and config
# ---------------------------------------------------------------------------
ROOT       <- "K:/ETH TOWERS"
GEE_DIR    <- file.path(ROOT, "ETH_towers")
TRAIN_CSV  <- file.path(ROOT, "training_table.csv")
RES_DIR    <- file.path(ROOT, "results")
dir_create(RES_DIR)

# Pipeline knobs you can edit
BLOCK_SIZE_M    <- 50000     # 50 km spatial blocks
N_OUTER_FOLDS   <- 10
N_INNER_FOLDS   <- 5
VIF_THRESHOLD   <- 10
COR_THRESHOLD   <- 0.85
TUNE_GRID_SIZE  <- 20        # candidates per model in tuning
SHAP_NSIM       <- 50        # SHAP Monte Carlo samples per row
DENSIFICATION_RADIUS_M <- 5000
SUITABILITY_CLASSES <- 5     # for Jenks classification

# Categorical predictors (treated as factors and not VIF-screened)
CATEGORICAL_NAMES <- c("ETH_landcover_2021", "ETH_smod_2020",
                       "ETH_protected_areas", "ETH_flood_binary")

# ---------------------------------------------------------------------------
# 2. Load training data
# ---------------------------------------------------------------------------
cat("=== 1. Load training table ===\n")
train_raw <- read_csv(TRAIN_CSV, show_col_types = FALSE)
cat("Loaded ", nrow(train_raw), " rows, ", ncol(train_raw), " cols\n", sep = "")
cat("Label distribution:\n"); print(table(train_raw$label))

# Identify predictor columns (everything except metadata)
META_COLS <- intersect(c("point_id", "label", "source",
                         "x_proj", "y_proj", "lon", "lat"),
                       names(train_raw))
predictor_cols <- setdiff(names(train_raw), META_COLS)
cat("Predictor columns: ", length(predictor_cols), "\n", sep = "")

# Drop any rows with NA in any predictor (rare; from script 02 these should be 0)
train <- train_raw %>% drop_na(all_of(predictor_cols))
cat("After NA-drop: ", nrow(train), " rows\n", sep = "")

# Coerce categorical predictors to factor where present
cat_cols  <- intersect(CATEGORICAL_NAMES, predictor_cols)
cont_cols <- setdiff(predictor_cols, cat_cols)
for (col in cat_cols) train[[col]] <- as.factor(train[[col]])

# Outcome must be a factor with "presence" as the positive (first) level
train$label_fct <- factor(ifelse(train$label == 1, "presence", "absence"),
                          levels = c("presence", "absence"))

# ---------------------------------------------------------------------------
# 3. Multicollinearity screening (VIF + Spearman) on continuous predictors
# ---------------------------------------------------------------------------
cat("\n=== 2. Multicollinearity screening ===\n")

screen_collinearity <- function(df, vars,
                                cor_thr = COR_THRESHOLD,
                                vif_thr = VIF_THRESHOLD) {
  # Iteratively drop the most collinear continuous predictor
  X <- df[, vars, drop = FALSE]
  X <- X[, sapply(X, is.numeric), drop = FALSE]
  dropped <- character(0)
  log_rows <- list()
  step <- 0
  repeat {
    step <- step + 1
    if (ncol(X) < 2) break
    # 1. Spearman correlation check
    cor_mat <- abs(cor(X, method = "spearman", use = "complete.obs"))
    diag(cor_mat) <- 0
    if (max(cor_mat, na.rm = TRUE) >= cor_thr) {
      worst_pair <- which(cor_mat == max(cor_mat), arr.ind = TRUE)[1, ]
      v1 <- colnames(X)[worst_pair[1]]; v2 <- colnames(X)[worst_pair[2]]
      mean_cor <- colMeans(cor_mat, na.rm = TRUE)
      drop <- if (mean_cor[v1] >= mean_cor[v2]) v1 else v2
      log_rows[[step]] <- tibble(step = step, rule = "rho",
                                 dropped = drop,
                                 value = round(cor_mat[v1, v2], 3))
      X[[drop]] <- NULL
      dropped <- c(dropped, drop)
      next
    }
    # 2. VIF check
    vifs <- sapply(seq_along(X), function(i) {
      y <- X[[i]]; Xo <- X[, -i, drop = FALSE]
      if (ncol(Xo) == 0) return(1)
      r2 <- summary(lm(y ~ ., data = Xo))$r.squared
      1 / (1 - r2)
    })
    names(vifs) <- colnames(X)
    if (max(vifs, na.rm = TRUE) >= vif_thr) {
      drop <- names(which.max(vifs))
      log_rows[[step]] <- tibble(step = step, rule = "vif",
                                 dropped = drop,
                                 value = round(max(vifs), 2))
      X[[drop]] <- NULL
      dropped <- c(dropped, drop)
    } else break
  }
  list(retained = colnames(X), dropped = dropped,
       log = bind_rows(log_rows))
}

screen <- screen_collinearity(train, cont_cols)
write_csv(screen$log, file.path(RES_DIR, "01_vif_screening.csv"))
retained_cont <- screen$retained
cat("Continuous predictors retained: ", length(retained_cont),
    " (dropped ", length(screen$dropped), ")\n", sep = "")
if (length(screen$dropped) > 0)
  cat("  Dropped: ", paste(screen$dropped, collapse = ", "), "\n", sep = "")
final_predictors <- c(retained_cont, cat_cols)
cat("Total predictors used: ", length(final_predictors), "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Spatial block cross-validation folds (blockCV)
# ---------------------------------------------------------------------------
cat("\n=== 3. Spatial block CV folds ===\n")

# Use projected coords already in the training table (x_proj, y_proj)
pts_sf <- st_as_sf(train, coords = c("x_proj", "y_proj"),
                   crs = 32637, remove = FALSE)

cv_blocks <- cv_spatial(
  x          = pts_sf,
  column     = "label",
  size       = BLOCK_SIZE_M,
  k          = N_OUTER_FOLDS,
  selection  = "random",
  iteration  = 50,
  biomod2    = FALSE,
  progress   = FALSE,
  plot       = FALSE,
  report     = FALSE,
  seed       = 42
)
train$fold <- cv_blocks$folds_ids
cat("Outer-fold sizes:\n"); print(table(train$fold))

# Build a rsample object so tidymodels can use these folds
fold_to_split <- function(k) {
  idx_test  <- which(train$fold == k)
  idx_train <- setdiff(seq_len(nrow(train)), idx_test)
  rsample::make_splits(
    list(analysis = idx_train, assessment = idx_test),
    data = train
  )
}
outer_resamples <- rsample::manual_rset(
  splits = lapply(seq_len(N_OUTER_FOLDS), fold_to_split),
  ids    = paste0("Fold", seq_len(N_OUTER_FOLDS))
)

# ---------------------------------------------------------------------------
# 5. Recipe (preprocessing pipeline)
# ---------------------------------------------------------------------------
cat("\n=== 4. Preprocessing recipe ===\n")

formula_str <- paste("label_fct ~", paste(final_predictors, collapse = " + "))
rec <- recipe(as.formula(formula_str), data = train) %>%
  step_novel(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

# ---------------------------------------------------------------------------
# 6. Model specifications
# ---------------------------------------------------------------------------
cat("\n=== 5. Model specifications ===\n")

# Random Forest (ranger)
spec_rf <- rand_forest(
  trees = 1000,
  mtry  = tune(),
  min_n = tune()
) %>%
  set_engine("ranger", importance = "permutation", num.threads = 4) %>%
  set_mode("classification")

# XGBoost
spec_xgb <- boost_tree(
  trees      = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  min_n      = tune(),
  sample_size = tune(),
  mtry       = tune(),
  loss_reduction = tune()
) %>%
  set_engine("xgboost", nthread = 4) %>%
  set_mode("classification")

# LightGBM (via bonsai)
spec_lgbm <- boost_tree(
  trees      = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  min_n      = tune(),
  mtry       = tune()
) %>%
  set_engine("lightgbm", num_threads = 4) %>%
  set_mode("classification")

# Workflows
wf_rf  <- workflow() %>% add_recipe(rec) %>% add_model(spec_rf)
wf_xgb <- workflow() %>% add_recipe(rec) %>% add_model(spec_xgb)
wf_lgbm<- workflow() %>% add_recipe(rec) %>% add_model(spec_lgbm)

# Tuning grids (Latin hypercube)
n_pred <- length(final_predictors)
mtry_max <- max(2, floor(sqrt(n_pred)))

set.seed(42)
grid_rf <- grid_latin_hypercube(
  finalize(mtry(range = c(2, mtry_max)), train[, final_predictors]),
  min_n(range = c(1, 20)),
  size = TUNE_GRID_SIZE
)
set.seed(43)
grid_xgb <- grid_latin_hypercube(
  trees(range = c(300, 1500)),
  tree_depth(range = c(3, 10)),
  learn_rate(range = c(-3, -1)),
  min_n(range = c(2, 25)),
  sample_prop(range = c(0.5, 1.0)),
  finalize(mtry(range = c(2, mtry_max)), train[, final_predictors]),
  loss_reduction(range = c(-3, 1)),
  size = TUNE_GRID_SIZE
)
set.seed(44)
grid_lgbm <- grid_latin_hypercube(
  trees(range = c(300, 1500)),
  tree_depth(range = c(3, 10)),
  learn_rate(range = c(-3, -1)),
  min_n(range = c(2, 25)),
  finalize(mtry(range = c(2, mtry_max)), train[, final_predictors]),
  size = TUNE_GRID_SIZE
)

# Metrics for tuning (we tune on AUC)
tune_metrics <- metric_set(roc_auc, mn_log_loss)

ctrl_grid_args <- control_stack_grid()  # required for stacking later

# ---------------------------------------------------------------------------
# 7. Tune each base learner on the spatial outer folds
# ---------------------------------------------------------------------------
cat("\n=== 6. Tune base learners (this is the slow step) ===\n")

cat("  Tuning Random Forest...\n")
res_rf <- tune_grid(wf_rf, resamples = outer_resamples, grid = grid_rf,
                    metrics = tune_metrics, control = ctrl_grid_args)

cat("  Tuning XGBoost...\n")
res_xgb <- tune_grid(wf_xgb, resamples = outer_resamples, grid = grid_xgb,
                     metrics = tune_metrics, control = ctrl_grid_args)

cat("  Tuning LightGBM...\n")
res_lgbm <- tune_grid(wf_lgbm, resamples = outer_resamples, grid = grid_lgbm,
                      metrics = tune_metrics, control = ctrl_grid_args)

# Best hyperparameters per learner (by AUC)
best_rf   <- select_best(res_rf,   metric = "roc_auc")
best_xgb  <- select_best(res_xgb,  metric = "roc_auc")
best_lgbm <- select_best(res_lgbm, metric = "roc_auc")

# Finalised workflows fit on full data (for prediction stage later)
wf_rf_final   <- finalize_workflow(wf_rf,   best_rf)
wf_xgb_final  <- finalize_workflow(wf_xgb,  best_xgb)
wf_lgbm_final <- finalize_workflow(wf_lgbm, best_lgbm)

# ---------------------------------------------------------------------------
# 8. Stacked ensemble (elastic-net meta-learner over RF + XGB + LGBM)
# ---------------------------------------------------------------------------
cat("\n=== 7. Stacked ensemble ===\n")
stk <- stacks() %>%
  add_candidates(res_rf,   name = "rf") %>%
  add_candidates(res_xgb,  name = "xgb") %>%
  add_candidates(res_lgbm, name = "lgbm") %>%
  blend_predictions(metric = metric_set(roc_auc),
                    penalty = 10 ^ seq(-4, -0.5, length.out = 20),
                    mixture = c(0.0, 0.25, 0.5, 0.75, 1.0)) %>%
  fit_members()

# ---------------------------------------------------------------------------
# 9. MaxEnt baseline (presence-background, fit per outer fold manually)
# ---------------------------------------------------------------------------
cat("\n=== 8. MaxEnt baseline ===\n")

if (!requireNamespace("maxnet", quietly = TRUE)) install.packages("maxnet")
library(maxnet)

# Build a numeric-only matrix for maxnet (it can't handle factors directly)
build_maxnet_matrix <- function(df) {
  X <- df[, final_predictors, drop = FALSE]
  for (col in cat_cols) {
    if (col %in% names(X))
      X[[col]] <- as.integer(as.factor(X[[col]]))
  }
  as.matrix(X)
}

# CV-fold predictions for MaxEnt
maxnet_cv_preds <- map_dfr(seq_len(N_OUTER_FOLDS), function(k) {
  idx_test  <- which(train$fold == k)
  idx_train <- setdiff(seq_len(nrow(train)), idx_test)
  X_tr <- build_maxnet_matrix(train[idx_train, ])
  y_tr <- as.integer(train$label[idx_train])
  X_te <- build_maxnet_matrix(train[idx_test, ])
  fit <- maxnet::maxnet(p = y_tr, data = as.data.frame(X_tr),
                        regmult = 1.0,
                        f = maxnet::maxnet.formula(p = y_tr,
                                                   data = as.data.frame(X_tr),
                                                   classes = "lqh"))
  preds <- predict(fit, newdata = as.data.frame(X_te), type = "logistic")
  tibble(model = "MaxEnt", fold = k,
         truth = factor(ifelse(train$label[idx_test] == 1,
                               "presence", "absence"),
                        levels = c("presence", "absence")),
         .pred_presence = as.numeric(preds))
})

# ---------------------------------------------------------------------------
# 10. Collect out-of-fold predictions for every model
# ---------------------------------------------------------------------------
cat("\n=== 9. Collect OOF predictions ===\n")

best_oof <- function(tune_res, best_pars, model_name) {
  collect_predictions(tune_res, parameters = best_pars,
                      summarize = FALSE) %>%
    transmute(model    = model_name,
              fold     = as.integer(gsub("Fold", "", id)),
              truth    = label_fct,
              .pred_presence = .pred_presence)
}

oof_rf   <- best_oof(res_rf,   best_rf,   "RandomForest")
oof_xgb  <- best_oof(res_xgb,  best_xgb,  "XGBoost")
oof_lgbm <- best_oof(res_lgbm, best_lgbm, "LightGBM")

# Stacked ensemble OOF: re-predict each fold's hold-out using the stack
oof_stack <- map_dfr(seq_len(N_OUTER_FOLDS), function(k) {
  idx_test <- which(train$fold == k)
  preds <- predict(stk, new_data = train[idx_test, ], type = "prob")
  tibble(model = "Stacked",
         fold  = k,
         truth = train$label_fct[idx_test],
         .pred_presence = preds$.pred_presence)
})

oof_all <- bind_rows(oof_rf, oof_xgb, oof_lgbm, oof_stack, maxnet_cv_preds)

# ---------------------------------------------------------------------------
# 11. Performance metrics per fold + summary across folds
# ---------------------------------------------------------------------------
cat("\n=== 10. Performance metrics ===\n")

# Helper: pick TSS-maximising threshold on a presence/absence vector
pick_threshold_tss <- function(truth_bin, p) {
  obs <- as.integer(truth_bin == "presence")
  thr_grid <- seq(0.05, 0.95, by = 0.01)
  best <- thr_grid[which.max(sapply(thr_grid, function(thr) {
    pred <- as.integer(p >= thr)
    tp <- sum(pred == 1 & obs == 1); fn <- sum(pred == 0 & obs == 1)
    tn <- sum(pred == 0 & obs == 0); fp <- sum(pred == 1 & obs == 0)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    sens + spec - 1
  }))]
  best
}

metrics_per_fold <- oof_all %>%
  group_by(model, fold) %>%
  group_modify(~ {
    df  <- .x
    obs <- as.integer(df$truth == "presence")
    p   <- df$.pred_presence
    if (length(unique(obs)) < 2) return(tibble())  # skip degenerate folds
    auc <- as.numeric(pROC::auc(pROC::roc(obs, p, quiet = TRUE)))
    auc_pr <- yardstick::pr_auc_vec(df$truth, p,
                                    event_level = "first")
    thr  <- pick_threshold_tss(df$truth, p)
    pred <- as.integer(p >= thr)
    tp <- sum(pred == 1 & obs == 1); fn <- sum(pred == 0 & obs == 1)
    tn <- sum(pred == 0 & obs == 0); fp <- sum(pred == 1 & obs == 0)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    tss  <- sens + spec - 1
    bal  <- (sens + spec) / 2
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
            2 * prec * sens / (prec + sens) else NA_real_
    n <- length(obs)
    po <- (tp + tn) / n
    pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
    kap <- if (pe < 1) (po - pe) / (1 - pe) else NA_real_
    tibble(auc = auc, auc_pr = auc_pr, tss = tss, kappa = kap,
           sens = sens, spec = spec, bal_acc = bal, f1 = f1, thr = thr)
  }) %>%
  ungroup()

write_csv(metrics_per_fold, file.path(RES_DIR, "02_cv_metrics_per_fold.csv"))

metrics_summary <- metrics_per_fold %>%
  pivot_longer(c(auc, auc_pr, tss, kappa, sens, spec, bal_acc, f1, thr),
               names_to = "metric", values_to = "value") %>%
  group_by(model, metric) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd   = sd(value,   na.rm = TRUE),
            n    = n(), .groups = "drop") %>%
  arrange(metric, desc(mean))

write_csv(metrics_summary, file.path(RES_DIR, "03_cv_metrics_summary.csv"))
cat("Saved per-fold + summary metrics.\n")
print(metrics_summary %>% filter(metric %in% c("auc", "tss", "f1")))

# ---------------------------------------------------------------------------
# 12. Wilcoxon paired tests: stacked ensemble vs each base learner
# ---------------------------------------------------------------------------
cat("\n=== 11. Wilcoxon paired tests ===\n")

wide_auc <- metrics_per_fold %>%
  select(model, fold, auc) %>%
  pivot_wider(names_from = model, values_from = auc)

wilc <- map_dfr(setdiff(unique(metrics_per_fold$model), "Stacked"),
                function(m) {
  if (!m %in% names(wide_auc)) return(NULL)
  test <- suppressWarnings(wilcox.test(wide_auc$Stacked, wide_auc[[m]],
                                       paired = TRUE))
  tibble(comparison = paste("Stacked vs", m),
         statistic  = unname(test$statistic),
         p_value    = test$p.value,
         mean_diff  = mean(wide_auc$Stacked - wide_auc[[m]], na.rm = TRUE))
})
write_csv(wilc, file.path(RES_DIR, "04_wilcoxon_pairwise.csv"))
print(wilc)

# ---------------------------------------------------------------------------
# 13. Variable importance for tree models
# ---------------------------------------------------------------------------
cat("\n=== 12. Variable importance ===\n")

fit_rf_full   <- fit(wf_rf_final,   data = train)
fit_xgb_full  <- fit(wf_xgb_final,  data = train)
fit_lgbm_full <- fit(wf_lgbm_final, data = train)

vi_rf <- vi(extract_fit_engine(fit_rf_full)) %>%
  mutate(model = "RandomForest")
vi_xgb <- vi(extract_fit_engine(fit_xgb_full)) %>%
  mutate(model = "XGBoost")
vi_lgbm <- vi(extract_fit_engine(fit_lgbm_full)) %>%
  mutate(model = "LightGBM")

write_csv(vi_rf,   file.path(RES_DIR, "05_var_importance_rf.csv"))
write_csv(vi_xgb,  file.path(RES_DIR, "05_var_importance_xgb.csv"))
write_csv(vi_lgbm, file.path(RES_DIR, "05_var_importance_lgbm.csv"))

# ---------------------------------------------------------------------------
# 14. SHAP values (TreeSHAP via fastshap on XGBoost)
# ---------------------------------------------------------------------------
cat("\n=== 13. SHAP (TreeSHAP via XGBoost) ===\n")

# Build the design matrix the way the XGBoost workflow saw it (post-recipe)
prep_rec <- prep(rec)
X_baked  <- bake(prep_rec, new_data = train, all_predictors())
X_mat    <- as.matrix(X_baked)

xgb_engine <- extract_fit_engine(fit_xgb_full)
shap_pred_fun <- function(object, newdata) {
  predict(object, newdata = as.matrix(newdata))
}
shap_vals <- fastshap::explain(xgb_engine, X = X_mat, nsim = SHAP_NSIM,
                               pred_wrapper = shap_pred_fun)
shap_df <- as_tibble(shap_vals)
shap_df$point_id <- if ("point_id" %in% names(train)) train$point_id else seq_len(nrow(train))
shap_df$label    <- train$label
write_csv(shap_df, file.path(RES_DIR, "06_shap_values.csv"))

# ---------------------------------------------------------------------------
# 15. National prediction surface (100 m)
# ---------------------------------------------------------------------------
cat("\n=== 14. National prediction surface ===\n")

# Load predictor stack (same files, same order, same names as training)
tif_files <- as.character(dir_ls(GEE_DIR, regexp = "\\.tif$", type = "file"))
tif_files <- tif_files[!grepl("^_", path_file(tif_files))]
stk_full  <- rast(tif_files)
names(stk_full) <- tools::file_path_sans_ext(path_file(tif_files))

# Subset stack to the predictors actually used by the model
need_layers <- final_predictors
missing_layers <- setdiff(need_layers, names(stk_full))
if (length(missing_layers) > 0) {
  warning("Missing layers in raster stack: ",
          paste(missing_layers, collapse = ", "))
}
stk_full <- stk_full[[intersect(need_layers, names(stk_full))]]

# Predict in chunks (100m national grids are huge)
predict_stack_with_workflow <- function(wf, raster_stack, type = "prob") {
  predict_fun <- function(model, ...) {
    df <- as.data.frame(...)
    # Make sure factor predictors get the same levels as in training
    for (col in cat_cols) {
      if (col %in% names(df)) {
        df[[col]] <- factor(df[[col]],
                            levels = levels(train[[col]]))
      }
    }
    p <- predict(model, new_data = df, type = type)
    p$.pred_presence
  }
  terra::predict(raster_stack, wf,
                 fun = predict_fun, na.rm = TRUE,
                 cores = 1, progress = "text")
}

# Use the stacked ensemble for the headline map
suit_continuous <- predict_stack_with_workflow(stk, stk_full, type = "prob")
names(suit_continuous) <- "suitability"
writeRaster(suit_continuous,
            file.path(RES_DIR, "07_suitability_continuous.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# Jenks classification into 5 classes
vals <- values(suit_continuous, na.rm = TRUE)
brks <- classIntervals(vals[seq(1, length(vals), length.out = 50000)],
                       n = SUITABILITY_CLASSES, style = "fisher")$brks
suit_classes <- classify(suit_continuous,
                         rcl = cbind(brks[-length(brks)], brks[-1],
                                     seq_len(SUITABILITY_CLASSES)),
                         include.lowest = TRUE, right = FALSE)
names(suit_classes) <- "suit_class"
writeRaster(suit_classes, file.path(RES_DIR, "07_suitability_classes.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# ---------------------------------------------------------------------------
# 16. Priority deployment surface (greenfield + densification)
# ---------------------------------------------------------------------------
cat("\n=== 15. Priority deployment surface ===\n")

rescale01 <- function(r) {
  v <- values(r, na.rm = TRUE)
  lo <- quantile(v, 0.01, na.rm = TRUE)
  hi <- quantile(v, 0.99, na.rm = TRUE)
  rr <- (r - lo) / (hi - lo)
  clamp(rr, 0, 1)
}

# Demand index = log-pop * built-up (both rescaled to 0-1)
pop_r <- rast(file.path(GEE_DIR, "ETH_population.tif"))
bld_r <- rast(file.path(GEE_DIR, "ETH_builtup_2020.tif"))
demand <- rescale01(log1p(pop_r)) * rescale01(bld_r)
demand <- resample(demand, suit_continuous, method = "bilinear")
demand <- rescale01(demand)

# Coverage gap = distance to existing tower, rescaled
dist_tower <- rast(file.path(GEE_DIR, "ETH_dist_to_existing_tower.tif"))
dist_tower <- resample(dist_tower, suit_continuous, method = "bilinear")
gap <- rescale01(dist_tower)

priority_score <- suit_continuous * demand * gap
names(priority_score) <- "priority"
writeRaster(priority_score, file.path(RES_DIR, "08_priority_score.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# Greenfield = top-decile priority AND >5km from existing tower
top_thr <- quantile(values(priority_score, na.rm = TRUE), 0.90, na.rm = TRUE)
greenfield <- (priority_score >= top_thr) &
              (dist_tower >= DENSIFICATION_RADIUS_M)
greenfield <- subst(greenfield, NA, 0)
names(greenfield) <- "priority_greenfield"
writeRaster(greenfield, file.path(RES_DIR, "08_priority_greenfield.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# Densification = high suitability + high demand + tower within 5km
densif <- (suit_continuous >= quantile(values(suit_continuous, na.rm = TRUE),
                                       0.75, na.rm = TRUE)) &
          (demand >= quantile(values(demand, na.rm = TRUE), 0.75,
                              na.rm = TRUE)) &
          (dist_tower < DENSIFICATION_RADIUS_M)
densif <- subst(densif, NA, 0)
names(densif) <- "priority_densification"
writeRaster(densif, file.path(RES_DIR, "08_priority_densification.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# ---------------------------------------------------------------------------
# 17. Run summary
# ---------------------------------------------------------------------------
cat("\n=== 16. Run summary ===\n")
sink(file.path(RES_DIR, "09_run_summary.txt"))
cat("Ethiopia Telecom Tower Siting - ML pipeline run summary\n")
cat("Generated: ", format(Sys.time()), "\n", sep = "")
cat("------------------------------------------------------\n\n")
cat("Training rows           : ", nrow(train), "\n", sep = "")
cat("Final predictors        : ", length(final_predictors), "\n", sep = "")
cat("  Continuous  : ", length(retained_cont), "\n", sep = "")
cat("  Categorical : ", length(cat_cols), "\n", sep = "")
cat("Predictors dropped (collinearity): ",
    paste(screen$dropped, collapse = ", "), "\n", sep = "")
cat("Spatial blocks: ", BLOCK_SIZE_M / 1000, " km, outer folds: ",
    N_OUTER_FOLDS, "\n", sep = "")
cat("\nBest hyperparameters:\n")
cat("  RF   :"); print(best_rf)
cat("  XGB  :"); print(best_xgb)
cat("  LGBM :"); print(best_lgbm)
cat("\nMetric summary (mean +/- sd across spatial folds):\n")
print(metrics_summary %>% filter(metric %in%
        c("auc", "auc_pr", "tss", "kappa", "f1", "bal_acc")), n = Inf)
cat("\nWilcoxon paired tests (stacked vs base learners, AUC):\n")
print(wilc)
sink()

cat("\n\nAll outputs in: ", RES_DIR, "\n", sep = "")
cat("Done.\n")
