# ============================================================================
# Ethiopia Telecom Tower Siting - Full ML Pipeline (v4 final)
# ----------------------------------------------------------------------------
# Critical fixes vs v3:
#   1. HONEST NESTED OOF STACKING.
#      Previous versions fit the stacked ensemble on the full training set
#      and then predicted on each held-out fold. Those predictions were not
#      truly out-of-sample (the stack had already seen every point during
#      its meta-learner fit), which is why the Stacked AUC came back at
#      0.9999 with sd 0.0002 - that's training-set performance, not
#      generalisation. v4 refits the entire stack inside each outer fold
#      using only that fold's training data, then predicts on the held-out
#      block. This is what Roberts et al. (2017) require for honest
#      reporting in spatial models.
#
#   2. POPULATION-WEIGHTED ABSENCES (handled in script 01c).
#      Pseudo-absences are now sampled with probability proportional to
#      log(population), forcing the classifier to discriminate between
#      similar populated locations rather than between cities and bushland.
#
# Inputs:
#   K:/ETH TOWERS/training_table.csv          (rebuild with script 02 v3
#                                              after running 01c)
#   K:/ETH TOWERS/ETH_towers/                 (predictor rasters)
#
# Outputs (K:/ETH TOWERS/results/ + figures/):
#   01_vif_screening.csv
#   02_cv_metrics_per_fold.csv
#   03_cv_metrics_summary.csv
#   04_wilcoxon_pairwise.csv
#   05_var_importance_*.csv
#   06_shap_values.csv
#   07_suitability_continuous.tif / classes.tif
#   08_priority_score.tif / greenfield.tif / densification.tif
#   09_run_summary.txt
#   figures/fig1..fig8.png  (all at 300 dpi)
#
# Runtime: 2-4 hours (longer than v3 because of nested OOF refits).
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Packages
# ---------------------------------------------------------------------------
required <- c(
  "terra", "sf", "dplyr", "tibble", "readr", "tidyr", "purrr", "fs",
  "tidymodels", "ranger", "xgboost", "lightgbm", "bonsai",
  "stacks", "yardstick",
  "blockCV", "fastshap", "vip", "pdp", "tidytext",
  "classInt", "pROC", "maxnet",
  "ggplot2", "patchwork", "tidyterra", "viridis", "scales"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(purrr); library(fs)
  library(tidymodels); library(ranger); library(xgboost)
  library(lightgbm); library(bonsai)
  library(stacks); library(yardstick)
  library(blockCV); library(fastshap); library(vip); library(pdp)
  library(tidytext)
  library(classInt); library(pROC); library(maxnet)
  library(ggplot2); library(patchwork); library(tidyterra)
  library(viridis); library(scales)
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
FIG_DIR    <- file.path(RES_DIR, "figures")
dir_create(RES_DIR); dir_create(FIG_DIR)

EXCLUDE_FROM_PREDICTORS <- c("ETH_ookla_mobile_dl_kbps",
                             "ETH_water_occurrence",
                             "ETH_dist_to_existing_tower")
EXCLUDE_FROM_RASTER_STACK <- c("ETH_ookla_mobile_dl_kbps",
                               "ETH_water_occurrence")
CATEGORICAL_NAMES <- c("ETH_landcover_2021", "ETH_smod_2020",
                       "ETH_protected_areas", "ETH_flood_binary")

BLOCK_SIZE_M    <- 50000
N_OUTER_FOLDS   <- 10
N_INNER_FOLDS   <- 5            # for stacking inside each outer fold
VIF_THRESHOLD   <- 10
COR_THRESHOLD   <- 0.85
TUNE_GRID_SIZE  <- 12           # smaller, since now tuned 10 x for nested OOF
SHAP_NSIM       <- 50
DENSIFICATION_RADIUS_M <- 5000
SUITABILITY_CLASSES <- 5

theme_pub <- function() {
  theme_bw(base_size = 11, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"))
}
save_fig <- function(p, file, w = 8, h = 5)
  ggsave(file, p, width = w, height = h, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# 2. Load training table
# ---------------------------------------------------------------------------
cat("=== 1. Load training table ===\n")
train_raw <- read_csv(TRAIN_CSV, show_col_types = FALSE)
cat("Rows: ", nrow(train_raw), "  Cols: ", ncol(train_raw), "\n", sep = "")
cat("Label distribution:\n"); print(table(train_raw$label))

META_COLS <- intersect(c("point_id", "label", "source",
                         "x_proj", "y_proj", "lon", "lat"),
                       names(train_raw))
predictor_cols <- setdiff(names(train_raw), META_COLS)
predictor_cols <- setdiff(predictor_cols, EXCLUDE_FROM_PREDICTORS)
cat("Predictor columns in use: ", length(predictor_cols), "\n", sep = "")
cat("Excluded from training: ",
    paste(EXCLUDE_FROM_PREDICTORS, collapse = ", "), "\n", sep = "")

train <- train_raw %>% drop_na(all_of(predictor_cols))
cat("After NA-drop: ", nrow(train), " rows\n", sep = "")

cat_cols  <- intersect(CATEGORICAL_NAMES, predictor_cols)
cont_cols <- setdiff(predictor_cols, cat_cols)
for (col in cat_cols) train[[col]] <- as.factor(train[[col]])

train$label_fct <- factor(ifelse(train$label == 1, "presence", "absence"),
                          levels = c("presence", "absence"))

# ---------------------------------------------------------------------------
# 3. Multicollinearity screening
# ---------------------------------------------------------------------------
cat("\n=== 2. Multicollinearity screening ===\n")
screen_collinearity <- function(df, vars, cor_thr, vif_thr) {
  X <- df[, vars, drop = FALSE]
  X <- X[, sapply(X, is.numeric), drop = FALSE]
  dropped <- character(0); log_rows <- list(); step <- 0
  repeat {
    step <- step + 1
    if (ncol(X) < 2) break
    cor_mat <- abs(cor(X, method = "spearman", use = "complete.obs"))
    diag(cor_mat) <- 0
    if (max(cor_mat, na.rm = TRUE) >= cor_thr) {
      worst <- which(cor_mat == max(cor_mat), arr.ind = TRUE)[1, ]
      v1 <- colnames(X)[worst[1]]; v2 <- colnames(X)[worst[2]]
      mc <- colMeans(cor_mat, na.rm = TRUE)
      drop <- if (mc[v1] >= mc[v2]) v1 else v2
      log_rows[[step]] <- tibble(step = step, rule = "rho", dropped = drop,
                                 value = round(cor_mat[v1, v2], 3))
      X[[drop]] <- NULL; dropped <- c(dropped, drop); next
    }
    vifs <- sapply(seq_along(X), function(i) {
      y <- X[[i]]; Xo <- X[, -i, drop = FALSE]
      if (ncol(Xo) == 0) return(1)
      r2 <- summary(lm(y ~ ., data = Xo))$r.squared
      1 / (1 - r2)
    })
    names(vifs) <- colnames(X)
    if (max(vifs, na.rm = TRUE) >= vif_thr) {
      drop <- names(which.max(vifs))
      log_rows[[step]] <- tibble(step = step, rule = "vif", dropped = drop,
                                 value = round(max(vifs), 2))
      X[[drop]] <- NULL; dropped <- c(dropped, drop)
    } else break
  }
  list(retained = colnames(X), dropped = dropped, log = bind_rows(log_rows))
}

screen <- screen_collinearity(train, cont_cols, COR_THRESHOLD, VIF_THRESHOLD)
write_csv(screen$log, file.path(RES_DIR, "01_vif_screening.csv"))
retained_cont    <- screen$retained
final_predictors <- c(retained_cont, cat_cols)
cat("Continuous retained: ", length(retained_cont),
    "  (dropped ", length(screen$dropped), ")\n", sep = "")
if (length(screen$dropped) > 0)
  cat("  Dropped: ", paste(screen$dropped, collapse = ", "), "\n", sep = "")
cat("Total predictors used: ", length(final_predictors), "\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Spatial block CV folds (outer)
# ---------------------------------------------------------------------------
cat("\n=== 3. Spatial block CV folds ===\n")
pts_sf <- st_as_sf(train, coords = c("x_proj", "y_proj"),
                   crs = 32637, remove = FALSE)
cv_blocks <- cv_spatial(
  x = pts_sf, column = "label",
  size = BLOCK_SIZE_M, k = N_OUTER_FOLDS,
  selection = "random", iteration = 50, biomod2 = FALSE,
  progress = FALSE, plot = FALSE, report = FALSE, seed = 42
)
train$fold <- cv_blocks$folds_ids
cat("Fold sizes:\n"); print(table(train$fold))

# Plot: spatial folds
p_folds <- ggplot() +
  geom_sf(data = st_transform(pts_sf, 4326),
          aes(colour = factor(cv_blocks$folds_ids)),
          size = 0.7, alpha = 0.85) +
  scale_colour_viridis_d(name = "Fold", option = "turbo") +
  labs(title = "Spatial block cross-validation folds",
       subtitle = paste0("Block size = ", BLOCK_SIZE_M / 1000, " km, ",
                         N_OUTER_FOLDS, " folds, ", nrow(train), " points")) +
  coord_sf() + theme_pub() +
  theme(panel.grid = element_line(colour = "grey95"))
save_fig(p_folds, file.path(FIG_DIR, "fig1_spatial_folds.png"), w = 8, h = 6)

# ---------------------------------------------------------------------------
# 5. Recipe and model specs (shared across nested OOF and final fits)
# ---------------------------------------------------------------------------
cat("\n=== 4. Recipe and model specs ===\n")
formula_str <- paste("label_fct ~", paste(final_predictors, collapse = " + "))
rec <- recipe(as.formula(formula_str), data = train) %>%
  step_novel(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

spec_rf <- rand_forest(trees = 1000, mtry = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation", num.threads = 4) %>%
  set_mode("classification")
spec_xgb <- boost_tree(trees = tune(), tree_depth = tune(),
                       learn_rate = tune(), min_n = tune(),
                       sample_size = tune(), mtry = tune(),
                       loss_reduction = tune()) %>%
  set_engine("xgboost", nthread = 4) %>%
  set_mode("classification")
spec_lgbm <- boost_tree(trees = tune(), tree_depth = tune(),
                        learn_rate = tune(), min_n = tune(),
                        mtry = tune()) %>%
  set_engine("lightgbm", num_threads = 4) %>%
  set_mode("classification")

wf_rf   <- workflow() %>% add_recipe(rec) %>% add_model(spec_rf)
wf_xgb  <- workflow() %>% add_recipe(rec) %>% add_model(spec_xgb)
wf_lgbm <- workflow() %>% add_recipe(rec) %>% add_model(spec_lgbm)

n_pred  <- length(final_predictors)
mtry_max <- max(2, floor(sqrt(n_pred)))

build_grids <- function(seed_offset = 0) {
  set.seed(42 + seed_offset)
  g_rf <- grid_latin_hypercube(
    finalize(mtry(range = c(2, mtry_max)), train[, final_predictors]),
    min_n(range = c(1, 20)), size = TUNE_GRID_SIZE)
  set.seed(43 + seed_offset)
  g_xgb <- grid_latin_hypercube(
    trees(range = c(300, 1500)), tree_depth(range = c(3, 10)),
    learn_rate(range = c(-3, -1)), min_n(range = c(2, 25)),
    sample_prop(range = c(0.5, 1.0)),
    finalize(mtry(range = c(2, mtry_max)), train[, final_predictors]),
    loss_reduction(range = c(-3, 1)), size = TUNE_GRID_SIZE)
  set.seed(44 + seed_offset)
  g_lgbm <- grid_latin_hypercube(
    trees(range = c(300, 1500)), tree_depth(range = c(3, 10)),
    learn_rate(range = c(-3, -1)), min_n(range = c(2, 25)),
    finalize(mtry(range = c(2, mtry_max)), train[, final_predictors]),
    size = TUNE_GRID_SIZE)
  list(rf = g_rf, xgb = g_xgb, lgbm = g_lgbm)
}

tune_metrics    <- metric_set(roc_auc, mn_log_loss)
ctrl_grid_args  <- control_stack_grid()

# ---------------------------------------------------------------------------
# 6. NESTED OUT-OF-FOLD STACKING (the honest computation)
# ---------------------------------------------------------------------------
cat("\n=== 5. Nested OOF stacking (slow, ~60-90 min) ===\n")
cat("    For each outer fold, refit a complete stack on the training\n")
cat("    blocks only, then predict on the held-out block. This is the\n")
cat("    only way to get honest spatial-CV metrics for the ensemble.\n\n")

run_one_outer_fold <- function(k) {
  cat(sprintf("--- Outer fold %d / %d ---\n", k, N_OUTER_FOLDS))
  idx_test  <- which(train$fold == k)
  idx_train <- setdiff(seq_len(nrow(train)), idx_test)
  train_k <- train[idx_train, ]
  test_k  <- train[idx_test, ]

  # Inner spatial folds within this outer training partition
  inner_pts <- st_as_sf(train_k, coords = c("x_proj", "y_proj"),
                        crs = 32637, remove = FALSE)
  inner_blocks <- tryCatch(
    cv_spatial(x = inner_pts, column = "label",
               size = BLOCK_SIZE_M, k = N_INNER_FOLDS,
               selection = "random", iteration = 30, biomod2 = FALSE,
               progress = FALSE, plot = FALSE, report = FALSE,
               seed = 100 + k),
    error = function(e) NULL
  )
  if (is.null(inner_blocks)) {
    # If blockCV fails for this fold (rare), fall back to random k-fold inner
    inner_ids <- sample(rep_len(seq_len(N_INNER_FOLDS), nrow(train_k)))
  } else {
    inner_ids <- inner_blocks$folds_ids
  }
  fold_to_split_inner <- function(j) {
    test_j  <- which(inner_ids == j)
    train_j <- setdiff(seq_len(nrow(train_k)), test_j)
    rsample::make_splits(list(analysis = train_j, assessment = test_j),
                         data = train_k)
  }
  inner_resamples <- rsample::manual_rset(
    splits = lapply(seq_len(N_INNER_FOLDS), fold_to_split_inner),
    ids    = paste0("InnerFold", seq_len(N_INNER_FOLDS))
  )

  grids <- build_grids(seed_offset = k)
  res_rf_k   <- tune_grid(wf_rf,   resamples = inner_resamples,
                          grid = grids$rf,
                          metrics = tune_metrics, control = ctrl_grid_args)
  res_xgb_k  <- tune_grid(wf_xgb,  resamples = inner_resamples,
                          grid = grids$xgb,
                          metrics = tune_metrics, control = ctrl_grid_args)
  res_lgbm_k <- tune_grid(wf_lgbm, resamples = inner_resamples,
                          grid = grids$lgbm,
                          metrics = tune_metrics, control = ctrl_grid_args)

  best_rf_k   <- select_best(res_rf_k,   metric = "roc_auc")
  best_xgb_k  <- select_best(res_xgb_k,  metric = "roc_auc")
  best_lgbm_k <- select_best(res_lgbm_k, metric = "roc_auc")

  # Base learners refit on full outer-train, predicted on held-out block
  fit_rf_k   <- fit(finalize_workflow(wf_rf,   best_rf_k),   data = train_k)
  fit_xgb_k  <- fit(finalize_workflow(wf_xgb,  best_xgb_k),  data = train_k)
  fit_lgbm_k <- fit(finalize_workflow(wf_lgbm, best_lgbm_k), data = train_k)

  pred_rf_k   <- predict(fit_rf_k,   test_k, type = "prob")$.pred_presence
  pred_xgb_k  <- predict(fit_xgb_k,  test_k, type = "prob")$.pred_presence
  pred_lgbm_k <- predict(fit_lgbm_k, test_k, type = "prob")$.pred_presence

  # Stacked ensemble fit ONLY on outer-train
  stack_k <- stacks() %>%
    add_candidates(res_rf_k,   name = "rf") %>%
    add_candidates(res_xgb_k,  name = "xgb") %>%
    add_candidates(res_lgbm_k, name = "lgbm") %>%
    blend_predictions(metric = metric_set(roc_auc),
                      penalty = 10 ^ seq(-4, -0.5, length.out = 15),
                      mixture = c(0.0, 0.25, 0.5, 0.75, 1.0)) %>%
    fit_members()
  pred_stack_k <- predict(stack_k, test_k, type = "prob")$.pred_presence

  # MaxEnt on outer-train
  build_maxnet_matrix <- function(df) {
    X <- df[, final_predictors, drop = FALSE]
    for (col in cat_cols) {
      if (col %in% names(X)) X[[col]] <- as.integer(as.factor(X[[col]]))
    }
    as.matrix(X)
  }
  X_tr <- build_maxnet_matrix(train_k); y_tr <- as.integer(train_k$label)
  X_te <- build_maxnet_matrix(test_k)
  fit_mx <- maxnet::maxnet(p = y_tr, data = as.data.frame(X_tr),
                           regmult = 1.0,
                           f = maxnet::maxnet.formula(p = y_tr,
                                  data = as.data.frame(X_tr),
                                  classes = "lqh"))
  pred_mx_k <- as.numeric(predict(fit_mx,
                                  newdata = as.data.frame(X_te),
                                  type = "logistic"))

  # Bind together
  truth_k <- factor(ifelse(train_k$label[idx_test - min(idx_test) + 1] == 1,
                           "presence", "absence"),
                    levels = c("presence", "absence"))
  # Above line is fragile across non-contiguous idx_test - rebuild from test_k
  truth_k <- test_k$label_fct

  bind_rows(
    tibble(model = "RandomForest", fold = k, truth = truth_k,
           .pred_presence = pred_rf_k),
    tibble(model = "XGBoost",      fold = k, truth = truth_k,
           .pred_presence = pred_xgb_k),
    tibble(model = "LightGBM",     fold = k, truth = truth_k,
           .pred_presence = pred_lgbm_k),
    tibble(model = "Stacked",      fold = k, truth = truth_k,
           .pred_presence = pred_stack_k),
    tibble(model = "MaxEnt",       fold = k, truth = truth_k,
           .pred_presence = pred_mx_k)
  )
}

oof_all <- map_dfr(seq_len(N_OUTER_FOLDS), run_one_outer_fold)
cat("Nested OOF predictions assembled: ", nrow(oof_all), " rows\n", sep = "")

# ---------------------------------------------------------------------------
# 7. Per-fold metrics + summary
# ---------------------------------------------------------------------------
cat("\n=== 6. Performance metrics ===\n")
pick_threshold_tss <- function(truth_bin, p) {
  obs <- as.integer(truth_bin == "presence")
  thr_grid <- seq(0.05, 0.95, by = 0.01)
  thr_grid[which.max(sapply(thr_grid, function(thr) {
    pred <- as.integer(p >= thr)
    tp <- sum(pred == 1 & obs == 1); fn <- sum(pred == 0 & obs == 1)
    tn <- sum(pred == 0 & obs == 0); fp <- sum(pred == 1 & obs == 0)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    sens + spec - 1
  }))]
}
metrics_per_fold <- oof_all %>%
  group_by(model, fold) %>%
  group_modify(~ {
    df <- .x
    obs <- as.integer(df$truth == "presence")
    p   <- df$.pred_presence
    if (length(unique(obs)) < 2) return(tibble())
    auc    <- as.numeric(pROC::auc(pROC::roc(obs, p, quiet = TRUE)))
    auc_pr <- yardstick::pr_auc_vec(df$truth, p, event_level = "first")
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
    n <- length(obs); po <- (tp + tn) / n
    pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
    kap <- if (pe < 1) (po - pe) / (1 - pe) else NA_real_
    tibble(auc = auc, auc_pr = auc_pr, tss = tss, kappa = kap,
           sens = sens, spec = spec, bal_acc = bal, f1 = f1, thr = thr)
  }) %>% ungroup()

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
print(metrics_summary %>% filter(metric %in% c("auc", "tss", "f1")))

# ROC curves and metric boxplots
roc_data <- oof_all %>%
  group_by(model) %>%
  group_modify(~ {
    obs <- as.integer(.x$truth == "presence")
    rc <- pROC::roc(obs, .x$.pred_presence, quiet = TRUE, direction = "<")
    tibble(fpr = 1 - rc$specificities, tpr = rc$sensitivities,
           auc = as.numeric(rc$auc))
  }) %>% ungroup()

auc_labels <- roc_data %>% group_by(model) %>%
  summarise(auc = first(auc), .groups = "drop") %>%
  mutate(label = sprintf("%s (AUC = %.3f)", model, auc))

p_roc <- ggplot(roc_data, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, linetype = "dashed", colour = "grey60") +
  geom_path(linewidth = 0.9) +
  scale_colour_viridis_d(name = NULL,
                         labels = setNames(auc_labels$label, auc_labels$model),
                         option = "D", end = 0.9) +
  labs(title = "ROC curves under spatial block cross-validation",
       subtitle = "Honest nested out-of-fold predictions",
       x = "False positive rate", y = "True positive rate") +
  coord_equal() + theme_pub()
save_fig(p_roc, file.path(FIG_DIR, "fig2_roc_curves.png"), w = 7, h = 6.3)

metric_long <- metrics_per_fold %>%
  pivot_longer(c(auc, tss, f1, kappa, bal_acc),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, auc = "AUC", tss = "TSS", f1 = "F1",
                         kappa = "Kappa", bal_acc = "Balanced acc."))
p_box <- ggplot(metric_long,
                aes(x = reorder(model, value, FUN = median),
                    y = value, fill = model)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.7, width = 0.6) +
  facet_wrap(~ metric, scales = "free_y", ncol = 5) +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(title = "Model performance across spatial CV folds", x = NULL, y = NULL) +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_box, file.path(FIG_DIR, "fig3_metric_boxplots.png"), w = 11, h = 4.5)

# ---------------------------------------------------------------------------
# 8. Wilcoxon paired tests
# ---------------------------------------------------------------------------
cat("\n=== 7. Wilcoxon paired tests ===\n")
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
# 9. Final fits on FULL training set (for prediction + interpretation)
# ---------------------------------------------------------------------------
cat("\n=== 8. Final fits on full training set ===\n")

outer_resamples <- rsample::manual_rset(
  splits = lapply(seq_len(N_OUTER_FOLDS), function(k) {
    test_k  <- which(train$fold == k)
    train_k <- setdiff(seq_len(nrow(train)), test_k)
    rsample::make_splits(list(analysis = train_k, assessment = test_k),
                         data = train)
  }),
  ids = paste0("Fold", seq_len(N_OUTER_FOLDS))
)

grids <- build_grids(seed_offset = 0)
cat("  Tuning RF (final)...\n")
res_rf_full   <- tune_grid(wf_rf,   resamples = outer_resamples,
                           grid = grids$rf,
                           metrics = tune_metrics, control = ctrl_grid_args)
cat("  Tuning XGB (final)...\n")
res_xgb_full  <- tune_grid(wf_xgb,  resamples = outer_resamples,
                           grid = grids$xgb,
                           metrics = tune_metrics, control = ctrl_grid_args)
cat("  Tuning LGBM (final)...\n")
res_lgbm_full <- tune_grid(wf_lgbm, resamples = outer_resamples,
                           grid = grids$lgbm,
                           metrics = tune_metrics, control = ctrl_grid_args)

best_rf   <- select_best(res_rf_full,   metric = "roc_auc")
best_xgb  <- select_best(res_xgb_full,  metric = "roc_auc")
best_lgbm <- select_best(res_lgbm_full, metric = "roc_auc")

fit_rf_full   <- fit(finalize_workflow(wf_rf,   best_rf),   data = train)
fit_xgb_full  <- fit(finalize_workflow(wf_xgb,  best_xgb),  data = train)
fit_lgbm_full <- fit(finalize_workflow(wf_lgbm, best_lgbm), data = train)

ens_stack <- stacks() %>%
  add_candidates(res_rf_full,   name = "rf") %>%
  add_candidates(res_xgb_full,  name = "xgb") %>%
  add_candidates(res_lgbm_full, name = "lgbm") %>%
  blend_predictions(metric = metric_set(roc_auc),
                    penalty = 10 ^ seq(-4, -0.5, length.out = 20),
                    mixture = c(0.0, 0.25, 0.5, 0.75, 1.0)) %>%
  fit_members()

# ---------------------------------------------------------------------------
# 10. Variable importance
# ---------------------------------------------------------------------------
cat("\n=== 9. Variable importance ===\n")
vi_rf   <- vi(extract_fit_engine(fit_rf_full))   %>% mutate(model = "RandomForest")
vi_xgb  <- vi(extract_fit_engine(fit_xgb_full))  %>% mutate(model = "XGBoost")
vi_lgbm <- vi(extract_fit_engine(fit_lgbm_full)) %>% mutate(model = "LightGBM")
write_csv(vi_rf,   file.path(RES_DIR, "05_var_importance_rf.csv"))
write_csv(vi_xgb,  file.path(RES_DIR, "05_var_importance_xgb.csv"))
write_csv(vi_lgbm, file.path(RES_DIR, "05_var_importance_lgbm.csv"))

vi_all <- bind_rows(vi_rf, vi_xgb, vi_lgbm) %>%
  group_by(model) %>%
  mutate(Importance = Importance / max(Importance)) %>%
  slice_max(order_by = Importance, n = 15) %>% ungroup()

p_vi <- ggplot(vi_all,
               aes(x = Importance,
                   y = reorder_within(Variable, Importance, model),
                   fill = model)) +
  geom_col(width = 0.7, alpha = 0.9) +
  facet_wrap(~ model, scales = "free_y") +
  scale_y_reordered() +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(title = "Top-15 predictors by permutation importance",
       x = "Importance (rescaled per model)", y = NULL) +
  theme_pub() + theme(legend.position = "none")
save_fig(p_vi, file.path(FIG_DIR, "fig4_var_importance.png"), w = 11, h = 6)

# ---------------------------------------------------------------------------
# 11. SHAP + plot
# ---------------------------------------------------------------------------
cat("\n=== 10. SHAP ===\n")
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

shap_long <- shap_df %>% select(-point_id, -label) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "shap")
shap_summary <- shap_long %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(shap), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_shap)) %>% slice_head(n = 15)

p_shap <- ggplot(shap_summary,
                 aes(x = mean_abs_shap, y = reorder(feature, mean_abs_shap))) +
  geom_col(fill = viridis(1, begin = 0.4), alpha = 0.9, width = 0.7) +
  labs(title = "Mean |SHAP value| - top 15 predictors (XGBoost)",
       x = "Mean |SHAP value|", y = NULL) +
  theme_pub()
save_fig(p_shap, file.path(FIG_DIR, "fig5_shap_summary.png"), w = 8, h = 6)

# ---------------------------------------------------------------------------
# 12. Partial dependence (top 4)
# ---------------------------------------------------------------------------
cat("\n=== 11. Partial dependence ===\n")
top_orig <- vi_xgb %>% arrange(desc(Importance)) %>%
  slice_head(n = 8) %>% pull(Variable) %>%
  intersect(final_predictors) %>% head(4)

pdp_data <- map_dfr(top_orig, function(v) {
  pd <- pdp::partial(fit_xgb_full, pred.var = v, train = train,
                     prob = TRUE, which.class = "presence",
                     grid.resolution = 30)
  tibble(feature = v, x = pd[[v]], y = pd[["yhat"]])
})

p_pdp <- ggplot(pdp_data, aes(x, y)) +
  geom_line(colour = viridis(1, begin = 0.4), linewidth = 0.9) +
  facet_wrap(~ feature, scales = "free_x", ncol = 2) +
  labs(title = "Partial dependence - top 4 predictors (XGBoost)",
       x = NULL, y = "Predicted P(tower)") +
  theme_pub()
save_fig(p_pdp, file.path(FIG_DIR, "fig6_partial_dependence.png"),
         w = 9, h = 6)

# ---------------------------------------------------------------------------
# 13. National prediction stack
# ---------------------------------------------------------------------------
cat("\n=== 12. National prediction stack ===\n")
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
pred_stack <- pred_stack_full[[intersect(final_predictors,
                                         names(pred_stack_full))]]
cat("Prediction stack layers: ", nlyr(pred_stack), "\n", sep = "")

predict_chunk <- function(model, ...) {
  df <- as.data.frame(...)
  for (col in cat_cols) {
    if (col %in% names(df))
      df[[col]] <- factor(df[[col]], levels = levels(train[[col]]))
  }
  preds <- predict(model, new_data = df, type = "prob")
  preds$.pred_presence
}

# ---------------------------------------------------------------------------
# 14. National 100m suitability surface
# ---------------------------------------------------------------------------
cat("\n=== 13. National suitability surface (20-40 min) ===\n")
suit_continuous <- terra::predict(pred_stack, ens_stack,
                                  fun = predict_chunk, na.rm = TRUE,
                                  cores = 1, progress = "text")
names(suit_continuous) <- "suitability"
writeRaster(suit_continuous,
            file.path(RES_DIR, "07_suitability_continuous.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

vals <- values(suit_continuous, na.rm = TRUE)
sample_idx <- seq(1, length(vals), length.out = min(50000, length(vals)))
brks <- classIntervals(vals[sample_idx], n = SUITABILITY_CLASSES,
                       style = "fisher")$brks
suit_classes <- classify(suit_continuous,
                         rcl = cbind(brks[-length(brks)], brks[-1],
                                     seq_len(SUITABILITY_CLASSES)),
                         include.lowest = TRUE, right = FALSE)
names(suit_classes) <- "suit_class"
writeRaster(suit_classes, file.path(RES_DIR, "07_suitability_classes.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

p_suit <- ggplot() +
  geom_spatraster(data = suit_continuous) +
  scale_fill_viridis_c(name = "Suitability\nP(tower)",
                       option = "magma", na.value = "transparent",
                       limits = c(0, 1)) +
  labs(title = "Predicted telecom tower suitability surface, Ethiopia",
       subtitle = "Stacked ensemble of RF + XGBoost + LightGBM") +
  coord_sf() + theme_pub() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 8))
save_fig(p_suit, file.path(FIG_DIR, "fig7_suitability_map.png"), w = 8, h = 7)

# ---------------------------------------------------------------------------
# 15. Priority deployment surface
# ---------------------------------------------------------------------------
cat("\n=== 14. Priority deployment surface ===\n")
rescale01 <- function(r) {
  v <- values(r, na.rm = TRUE)
  lo <- quantile(v, 0.01, na.rm = TRUE)
  hi <- quantile(v, 0.99, na.rm = TRUE)
  rr <- (r - lo) / (hi - lo)
  clamp(rr, 0, 1)
}

pop_r <- rast(file.path(GEE_DIR, "ETH_population.tif"))
bld_r <- rast(file.path(GEE_DIR, "ETH_builtup_2020.tif"))
demand <- rescale01(log1p(pop_r)) * rescale01(bld_r)
demand <- resample(demand, suit_continuous, method = "bilinear")
demand <- rescale01(demand)

dist_tower <- rast(file.path(GEE_DIR, "ETH_dist_to_existing_tower.tif"))
dist_tower <- resample(dist_tower, suit_continuous, method = "bilinear")
gap <- rescale01(dist_tower)

priority_score <- suit_continuous * demand * gap
names(priority_score) <- "priority"
writeRaster(priority_score, file.path(RES_DIR, "08_priority_score.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

top_thr <- quantile(values(priority_score, na.rm = TRUE), 0.90, na.rm = TRUE)
greenfield <- (priority_score >= top_thr) &
              (dist_tower >= DENSIFICATION_RADIUS_M)
greenfield <- subst(greenfield, NA, 0)
names(greenfield) <- "priority_greenfield"
writeRaster(greenfield, file.path(RES_DIR, "08_priority_greenfield.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

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

priority_cat <- ifel(greenfield == 1, 2,
                     ifel(densif == 1, 1, 0))
names(priority_cat) <- "priority_class"

p_pri <- ggplot() +
  geom_spatraster(data = priority_cat) +
  scale_fill_manual(name = "Priority class",
                    values = c("0" = "grey90", "1" = "#E07B00",
                               "2" = "#1B7837"),
                    labels = c("0" = "Other", "1" = "Densification",
                               "2" = "Greenfield"),
                    na.value = "transparent", drop = FALSE) +
  labs(title = "Priority deployment surface, Ethiopia",
       subtitle = "Greenfield = top-decile priority + >5 km from existing tower\nDensification = high suitability + high demand + within 5 km") +
  coord_sf() + theme_pub() +
  theme(panel.grid = element_blank())
save_fig(p_pri, file.path(FIG_DIR, "fig8_priority_map.png"), w = 8, h = 7)

# ---------------------------------------------------------------------------
# 16. Final run summary
# ---------------------------------------------------------------------------
cat("\n=== 15. Run summary ===\n")
sink(file.path(RES_DIR, "09_run_summary.txt"))
cat("Ethiopia Telecom Tower Siting - ML pipeline run summary (v4)\n")
cat("Generated: ", format(Sys.time()), "\n", sep = "")
cat("------------------------------------------------------------\n\n")
cat("Training rows           : ", nrow(train), "\n", sep = "")
cat("Final predictors        : ", length(final_predictors), "\n", sep = "")
cat("  Continuous  : ", length(retained_cont), "\n", sep = "")
cat("  Categorical : ", length(cat_cols), "\n", sep = "")
cat("Predictors excluded ex ante (leakage / no signal):\n  ",
    paste(EXCLUDE_FROM_PREDICTORS, collapse = ", "), "\n", sep = "")
cat("Predictors dropped by VIF / |rho|: ",
    paste(screen$dropped, collapse = ", "), "\n", sep = "")
cat("Spatial blocks: ", BLOCK_SIZE_M / 1000, " km, outer folds: ",
    N_OUTER_FOLDS, ", inner folds: ", N_INNER_FOLDS, "\n", sep = "")
cat("\nBest hyperparameters (final fits):\n")
cat("  RF   :"); print(best_rf)
cat("  XGB  :"); print(best_xgb)
cat("  LGBM :"); print(best_lgbm)
cat("\nMetric summary - HONEST nested OOF (mean +/- sd over folds):\n")
print(metrics_summary %>% filter(metric %in%
        c("auc", "auc_pr", "tss", "kappa", "f1", "bal_acc")), n = Inf)
cat("\nWilcoxon paired tests (stacked vs base learners, AUC):\n")
print(wilc)
sink()

cat("\nAll outputs in: ", RES_DIR, "\n", sep = "")
cat("All figures in: ", FIG_DIR, "\n", sep = "")
cat("DONE.\n")
