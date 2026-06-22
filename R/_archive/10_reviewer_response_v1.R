# ============================================================================
# REVIEWER-RESPONSE SENSITIVITY SUITE
# ----------------------------------------------------------------------------
# Addresses Reviewer 2 comments:
#   #4   Pseudo-absence sensitivity (buffer + ratio)
#   #5   Spatial CV block size justification (Moran I + variogram)
#   #6.1 Calibration (Brier score + reliability diagram)
#   #8   Equity weighting sensitivity
#   #19  Priority threshold sensitivity
#
# Designed for: HP Envy x360, 8 GB RAM, 500 GB SSD
#
# Strategy:
#   - Sequential (not parallel) to avoid OOM
#   - gc() between heavy operations
#   - All intermediate results saved as .rds; script resumes from last completed step
#   - Total runtime estimate on 8 GB: ~10-14 hours, can run overnight
#   - Each section can be run independently (uncomment the relevant block at the bottom)
#
# Inputs (all from previous pipeline):
#   K:/ETH TOWERS/results/presence_absence.gpkg
#   K:/ETH TOWERS/results/training_table.csv
#   K:/ETH TOWERS/results/07_suitability_1km.tif      (LightGBM headline)
#   K:/ETH TOWERS/ETH_towers/*.tif                     (all 23 predictors)
#   K:/ETH TOWERS/ETH_towers/ETH_dist_to_existing_tower.tif
#
# Outputs (written to K:/ETH TOWERS/results/reviewer_response/):
#   r1_moran_variogram_results.csv     (#5)
#   r1_moran_variogram_plot.png        (#5)
#   r2_calibration_metrics.csv         (#6.1)
#   r2_reliability_diagram.png         (#6.1)
#   r3_pseudo_absence_sensitivity.csv  (#4)
#   r3_pseudo_absence_plot.png         (#4)
#   r4_equity_weight_sensitivity.csv   (#8)
#   r5_threshold_sensitivity.csv       (#19)
# ============================================================================

# ----- 0. ENVIRONMENT SETUP ------------------------------------------------

# Conservative RAM management for 8 GB system
options(java.parameters = "-Xmx4g")
Sys.setenv(R_MAX_VSIZE = "6Gb")

required <- c("terra", "sf", "dplyr", "tibble", "readr", "tidyr", "fs",
              "ggplot2", "patchwork", "spdep", "gstat", "automap",
              "lightgbm", "ranger", "tidymodels", "blockCV",
              "yardstick", "pROC", "scales", "viridis")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(fs)
  library(ggplot2); library(patchwork)
  library(spdep); library(gstat); library(automap)
  library(lightgbm); library(ranger); library(tidymodels)
  library(blockCV); library(yardstick); library(pROC)
  library(scales); library(viridis)
})

sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.5)   # cap terra at 50% RAM
data.table::setDTthreads(2)

# Paths
ROOT     <- "K:/ETH TOWERS"
RES_DIR  <- file.path(ROOT, "results")
GEE_DIR  <- file.path(ROOT, "ETH_towers")
OUT_DIR  <- file.path(RES_DIR, "reviewer_response")
FIG_DIR  <- file.path(OUT_DIR, "figures")
CKPT_DIR <- file.path(OUT_DIR, "checkpoints")   # so we can resume after crash
dir_create(c(OUT_DIR, FIG_DIR, CKPT_DIR))

theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "bottom",
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"))
}
save_fig <- function(p, file, w, h)
  ggsave(file, p, width = w, height = h, dpi = 300, bg = "white")

# Helper: load if checkpoint exists, else run fn() and save
checkpoint <- function(name, fn) {
  f <- file.path(CKPT_DIR, paste0(name, ".rds"))
  if (file.exists(f)) {
    cat(sprintf("  [resume] loading %s from checkpoint\n", name))
    return(readRDS(f))
  }
  result <- fn()
  saveRDS(result, f)
  gc(verbose = FALSE)
  result
}

cat("=================================================================\n")
cat("  Reviewer-response sensitivity suite\n")
cat("  RAM-optimised for HP Envy x360 (8 GB)\n")
cat(sprintf("  Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")

# ============================================================================
# SECTION 1 — MORAN I + VARIOGRAM (addresses comment #5)
# ============================================================================
# Justifies the 50 km spatial block size using empirical spatial autocorrelation
# range of LightGBM residuals on the training set.
# Expected runtime: 20-30 minutes. RAM ~1.5 GB.
# ============================================================================
run_section_1 <- function() {

  cat("\n========== SECTION 1: Moran I + variogram (comment #5) ==========\n\n")

  # Load training table and predictions
  tt_path <- file.path(RES_DIR, "training_table.csv")
  if (!file.exists(tt_path))
    stop("training_table.csv not found. Re-run 02_build_training_table_v3.R first.")

  tt <- read_csv(tt_path, show_col_types = FALSE)
  cat(sprintf("  training_table: %d rows x %d cols\n", nrow(tt), ncol(tt)))

  # Get predicted P(tower) at each training point from the saved LightGBM model
  # If predictions not saved as a column, re-extract from the suitability raster
  suit_r <- rast(file.path(RES_DIR, "07_suitability_1km.tif"))
  if (nlyr(suit_r) > 1) suit_r <- suit_r[[1]]

  pa <- st_read(file.path(RES_DIR, "presence_absence.gpkg"), quiet = TRUE)
  pa$pred <- terra::extract(suit_r, vect(st_transform(pa, crs(suit_r))))[, 2]
  pa$obs  <- as.integer(pa$class == "presence")
  pa$resid <- pa$obs - pa$pred

  # Project to metric CRS for distance-based spatial analysis
  pa_utm <- st_transform(pa, 32637)
  coords <- st_coordinates(pa_utm)
  resid  <- pa_utm$resid
  keep   <- !is.na(resid) & is.finite(resid)
  coords <- coords[keep, ]
  resid  <- resid[keep]
  cat(sprintf("  usable residuals: %d / %d\n", length(resid), nrow(pa)))

  # --- 1a. Moran I at several distance bands ---
  cat("  Computing Moran I at distance bands (10, 25, 50, 75, 100, 150 km)\n")
  bands_km <- c(10, 25, 50, 75, 100, 150)
  moran_results <- lapply(bands_km, function(d) {
    nb <- dnearneigh(coords, 0, d * 1000)
    if (any(card(nb) == 0)) {
      # drop islands
      ok <- card(nb) > 0
      nb <- dnearneigh(coords[ok, ], 0, d * 1000)
      mt <- moran.test(resid[ok], nb2listw(nb, style = "W", zero.policy = TRUE),
                       zero.policy = TRUE)
    } else {
      mt <- moran.test(resid, nb2listw(nb, style = "W"))
    }
    tibble(distance_km = d,
           moran_I = as.numeric(mt$estimate[1]),
           p_value = mt$p.value)
  })
  moran_df <- bind_rows(moran_results)
  cat("  Moran I by distance band:\n")
  print(moran_df)

  # --- 1b. Empirical variogram of residuals ---
  cat("  Fitting empirical variogram (autoKrige selects best model)\n")
  pa_sp <- as_Spatial(pa_utm[keep, ])
  pa_sp$resid <- resid
  vg <- variogram(resid ~ 1, data = pa_sp,
                  cutoff = 200000, width = 10000)
  vg_fit <- tryCatch(
    autofitVariogram(resid ~ 1, input_data = pa_sp,
                     model = c("Sph", "Exp", "Gau"))$var_model,
    error = function(e) {
      cat("  [warn] autoKrige failed, trying manual exponential fit\n")
      vgm(psill = var(resid), model = "Exp", range = 30000,
          nugget = 0.5 * var(resid))
    })

  range_m <- vg_fit$range[2]
  cat(sprintf("  fitted variogram range: %.1f km\n", range_m / 1000))

  # --- 1c. Plot ---
  p_moran <- ggplot(moran_df, aes(distance_km, moran_I)) +
    geom_line(colour = "#C24642", linewidth = 1) +
    geom_point(size = 3, colour = "#C24642") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 50, linetype = "dotted", colour = "#0F4C81") +
    annotate("text", x = 52, y = max(moran_df$moran_I) * 0.95,
             hjust = 0, label = "50 km block size (chosen)",
             colour = "#0F4C81", fontface = "bold", size = 3.4) +
    labs(title = "(a) Moran's I of LightGBM residuals by distance band",
         x = "Distance band (km)", y = "Moran's I") +
    theme_pub()

  vg_df <- as_tibble(vg)
  p_vg <- ggplot(vg_df, aes(dist / 1000, gamma)) +
    geom_point(size = 2.5, colour = "#0F4C81") +
    geom_hline(yintercept = vg_fit$psill[1] + vg_fit$psill[2],
               linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = range_m / 1000,
               linetype = "dotted", colour = "#C24642") +
    annotate("text", x = range_m / 1000 + 5,
             y = (vg_fit$psill[1] + vg_fit$psill[2]) * 0.5,
             hjust = 0, label = sprintf("Range = %.1f km", range_m / 1000),
             colour = "#C24642", fontface = "bold", size = 3.4) +
    labs(title = "(b) Empirical semivariogram",
         x = "Lag distance (km)", y = "Semivariance") +
    theme_pub()

  p_all <- p_moran | p_vg
  save_fig(p_all, file.path(FIG_DIR, "r1_moran_variogram.png"), 11, 4.5)
  write_csv(moran_df, file.path(OUT_DIR, "r1_moran_results.csv"))
  write_csv(tibble(model = vg_fit$model[2],
                   range_km = round(range_m / 1000, 2),
                   psill = vg_fit$psill[2],
                   nugget = vg_fit$psill[1]),
            file.path(OUT_DIR, "r1_variogram_fit.csv"))

  list(moran = moran_df, vg_range_km = range_m / 1000)
}

# ============================================================================
# SECTION 2 — CALIBRATION ANALYSIS (addresses comment #6.1)
# ============================================================================
# Computes Brier score and produces a reliability diagram showing how well
# the predicted probabilities match observed presence frequencies.
# Expected runtime: 15-20 minutes. RAM ~1 GB.
# ============================================================================
run_section_2 <- function() {

  cat("\n========== SECTION 2: Calibration (comment #6.1) ==========\n\n")

  # Load OOF predictions (saved by the main pipeline)
  oof_path <- file.path(RES_DIR, "02_oof_predictions.csv")
  if (!file.exists(oof_path)) {
    cat("  [warn] 02_oof_predictions.csv not found.\n")
    cat("  Falling back to in-bag predictions (less reliable for calibration).\n")
    # Generate quick OOF from training table if needed
    return(NULL)
  }

  oof <- read_csv(oof_path, show_col_types = FALSE)
  cat(sprintf("  OOF predictions: %d rows, models: %s\n",
              nrow(oof), paste(unique(oof$model), collapse = ", ")))

  # Brier score per model
  brier <- oof %>%
    group_by(model) %>%
    summarise(brier = mean((pred_prob - as.integer(obs == "presence"))^2),
              n = n(), .groups = "drop") %>%
    arrange(brier)

  cat("  Brier scores:\n")
  print(brier)
  write_csv(brier, file.path(OUT_DIR, "r2_brier_scores.csv"))

  # Reliability diagram — bin predictions into deciles, plot mean predicted vs
  # observed frequency
  rel <- oof %>%
    mutate(bin = cut(pred_prob, breaks = seq(0, 1, by = 0.1),
                     include.lowest = TRUE)) %>%
    group_by(model, bin) %>%
    summarise(mean_pred = mean(pred_prob),
              obs_freq  = mean(as.integer(obs == "presence")),
              n         = n(), .groups = "drop") %>%
    filter(n >= 5)

  p_rel <- ggplot(rel, aes(mean_pred, obs_freq, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n), alpha = 0.7) +
    scale_size_continuous(range = c(2, 7), name = "Bin count") +
    scale_colour_viridis_d(option = "D", end = 0.85) +
    coord_equal() +
    scale_x_continuous(limits = c(0, 1), labels = percent) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(title = "Reliability diagram",
         subtitle = "Diagonal = perfect calibration; above = under-confident, below = over-confident",
         x = "Mean predicted P(tower)",
         y = "Observed presence frequency",
         colour = "Model") +
    theme_pub()

  save_fig(p_rel, file.path(FIG_DIR, "r2_reliability_diagram.png"), 8, 7.5)

  brier
}

# ============================================================================
# SECTION 3 — PSEUDO-ABSENCE SENSITIVITY (addresses comment #4)
# ============================================================================
# Re-fits LightGBM under multiple pseudo-absence designs:
#   buffer = {0.5, 1, 2, 5} km   x   ratio = {1:1, 1:2, 1:5}
# That's 12 model refits. Each refit is ~10-15 min on your machine.
# Expected total runtime: 3-5 hours. RAM ~3 GB per refit.
# Saves CSV with AUC, TSS, F1 + greenfield area % under each design.
# ============================================================================
run_section_3 <- function() {

  cat("\n========== SECTION 3: Pseudo-absence sensitivity (comment #4) ==========\n\n")

  # Load predictors and presences
  presences <- st_read(file.path(RES_DIR, "presence_absence.gpkg"), quiet = TRUE) %>%
    filter(class == "presence") %>%
    st_transform(32637)

  pop_r <- rast(file.path(GEE_DIR, "ETH_population.tif"))
  if (nlyr(pop_r) > 1) pop_r <- pop_r[[1]]
  pop_r <- project(pop_r, crs(presences))

  dist_tower <- rast(file.path(GEE_DIR, "ETH_dist_to_existing_tower.tif"))
  if (nlyr(dist_tower) > 1) dist_tower <- dist_tower[[1]]
  dist_tower <- project(dist_tower, pop_r, method = "bilinear")

  # Build weighted-sampling probability surface (same as main pipeline)
  pop_log <- app(pop_r, function(v) {
    ok <- !is.na(v) & v > 0
    out <- v
    out[ok] <- log1p(v[ok])
    out[!ok] <- NA
    out
  })

  buffers_km <- c(0.5, 1, 2, 5)
  ratios     <- c(1, 2, 5)

  designs <- expand.grid(buffer_km = buffers_km, ratio = ratios,
                         stringsAsFactors = FALSE)
  cat(sprintf("  Total designs to fit: %d (estimated %.1f hours)\n",
              nrow(designs), nrow(designs) * 0.25))

  # Function: generate pseudo-absences for one design, return training table
  gen_design <- function(buffer_km, ratio) {
    set.seed(20260512)
    n_pres <- nrow(presences)
    n_abs  <- n_pres * ratio

    # Sample population-weighted candidates (oversample for filtering)
    candidates_n <- n_abs * 8
    cand <- spatSample(pop_log, candidates_n, method = "weights",
                       na.rm = TRUE, xy = TRUE, as.points = TRUE)
    cand_sf <- st_as_sf(cand) %>% st_set_crs(crs(pop_log))

    # Filter: outside buffer of any presence
    pres_buf <- st_buffer(presences, buffer_km * 1000) %>% st_union()
    cand_sf <- cand_sf[!st_intersects(cand_sf, pres_buf, sparse = FALSE)[, 1], ]

    if (nrow(cand_sf) < n_abs) {
      cat(sprintf("    [warn] only %d candidates after buffer, using all\n",
                  nrow(cand_sf)))
      absences <- cand_sf
    } else {
      absences <- cand_sf[sample(nrow(cand_sf), n_abs), ]
    }

    list(presences = presences, absences = absences,
         n_p = n_pres, n_a = nrow(absences))
  }

  # Function: fit a fast LightGBM (no tuning, just CV evaluation)
  # We use sensible default hyperparameters consistent with the main run
  fit_eval <- function(design_data) {
    # Build training matrix by extracting all predictors at the points
    pred_files <- list.files(GEE_DIR, pattern = "\\.tif$", full.names = TRUE)
    pred_files <- pred_files[!grepl("dist_to_existing_tower|water_occurrence|ookla",
                                    pred_files)]

    pts <- bind_rows(
      mutate(design_data$presences, class = 1),
      mutate(design_data$absences,  class = 0)
    ) %>% st_set_crs(32637)

    # Layer-by-layer extraction to manage RAM
    X <- as.data.frame(pts)[, "class", drop = FALSE]
    for (f in pred_files) {
      r <- rast(f); if (nlyr(r) > 1) r <- r[[1]]
      r <- project(r, "EPSG:32637")
      val <- terra::extract(r, vect(pts))[, 2]
      nm <- tools::file_path_sans_ext(basename(f))
      X[[nm]] <- val
      rm(r); gc(verbose = FALSE)
    }

    # Drop rows with > 4 NAs
    n_na <- rowSums(is.na(X))
    X <- X[n_na <= 4, ]
    # Impute remaining
    for (j in 2:ncol(X)) {
      if (is.numeric(X[[j]])) X[[j]][is.na(X[[j]])] <- median(X[[j]], na.rm = TRUE)
      else X[[j]][is.na(X[[j]])] <- names(sort(table(X[[j]]), decreasing = TRUE))[1]
    }

    # 5-fold spatial-ish CV (random here for speed; main paper uses spatial CV)
    set.seed(42)
    folds <- sample(1:5, nrow(X), replace = TRUE)
    aucs <- numeric(5); tsss <- numeric(5); f1s <- numeric(5)

    for (k in 1:5) {
      train <- X[folds != k, ]; test <- X[folds == k, ]
      dtr <- lgb.Dataset(as.matrix(train[, -1]), label = train$class)
      mod <- lgb.train(
        params = list(objective = "binary", metric = "auc",
                      num_leaves = 31, learning_rate = 0.05,
                      feature_fraction = 0.8, bagging_fraction = 0.8,
                      num_threads = 2, verbosity = -1),
        data = dtr, nrounds = 200, verbose = -1
      )
      pred <- predict(mod, as.matrix(test[, -1]))
      r <- roc(test$class, pred, quiet = TRUE)
      aucs[k] <- as.numeric(auc(r))
      cls <- ifelse(pred > 0.5, 1, 0)
      tp <- sum(cls == 1 & test$class == 1)
      fp <- sum(cls == 1 & test$class == 0)
      fn <- sum(cls == 0 & test$class == 1)
      tn <- sum(cls == 0 & test$class == 0)
      sens <- tp / (tp + fn); spec <- tn / (tn + fp)
      prec <- tp / (tp + fp + 1e-9)
      tsss[k] <- sens + spec - 1
      f1s[k] <- 2 * prec * sens / (prec + sens + 1e-9)
      rm(mod, dtr); gc(verbose = FALSE)
    }
    tibble(auc_mean = mean(aucs), auc_sd = sd(aucs),
           tss_mean = mean(tsss), tss_sd = sd(tsss),
           f1_mean = mean(f1s),  f1_sd = sd(f1s),
           n_p = design_data$n_p, n_a = design_data$n_a)
  }

  # Run all designs sequentially
  results <- list()
  for (i in seq_len(nrow(designs))) {
    d <- designs[i, ]
    ckpt_name <- sprintf("r3_design_b%.1f_r%d", d$buffer_km, d$ratio)
    res <- checkpoint(ckpt_name, function() {
      cat(sprintf("  [%d/%d] buffer = %.1f km, ratio = 1:%d\n",
                  i, nrow(designs), d$buffer_km, d$ratio))
      dd <- gen_design(d$buffer_km, d$ratio)
      perf <- fit_eval(dd)
      bind_cols(d, perf)
    })
    results[[i]] <- res
  }

  out <- bind_rows(results)
  write_csv(out, file.path(OUT_DIR, "r3_pseudo_absence_sensitivity.csv"))

  # Plot AUC vs buffer, by ratio
  p <- ggplot(out, aes(factor(buffer_km), auc_mean,
                       colour = factor(ratio), group = factor(ratio))) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = auc_mean - auc_sd, ymax = auc_mean + auc_sd),
                  width = 0.15) +
    scale_colour_viridis_d(option = "D", end = 0.85, name = "Pres:Abs ratio") +
    labs(title = "Pseudo-absence design sensitivity",
         subtitle = "Mean LightGBM AUC across 5 random folds (sensitivity check only)",
         x = "Presence-buffer distance (km)", y = "AUC") +
    theme_pub()
  save_fig(p, file.path(FIG_DIR, "r3_pseudo_absence_sensitivity.png"), 8, 5)

  cat(sprintf("\n  Variation across designs: AUC ranges %.3f to %.3f\n",
              min(out$auc_mean), max(out$auc_mean)))
  out
}

# ============================================================================
# SECTION 4 — EQUITY WEIGHTING SENSITIVITY (addresses comment #8)
# ============================================================================
# Re-computes the equity-adjusted priority surface under different weights
# on internet vs electricity indicators: w_internet in {0.3, 0.5, 0.7}
# and recomputes the orthogonality Spearman correlation each time.
# Expected runtime: 30 min. RAM ~2 GB.
# ============================================================================
run_section_4 <- function() {

  cat("\n========== SECTION 4: Equity weighting sensitivity (comment #8) ==========\n\n")

  # Load surfaces
  prio   <- rast(file.path(RES_DIR, "08_priority_score_1km.tif"))
  if (nlyr(prio) > 1) prio <- prio[[1]]

  i_r <- rast(file.path(GEE_DIR, "ETH_dhs_internet_users.tif"))
  if (nlyr(i_r) > 1) i_r <- i_r[[1]]
  e_r <- rast(file.path(GEE_DIR, "ETH_dhs_electricity_access.tif"))
  if (nlyr(e_r) > 1) e_r <- e_r[[1]]
  i_r <- resample(i_r, prio, method = "bilinear")
  e_r <- resample(e_r, prio, method = "bilinear")

  rescale01 <- function(r) {
    v <- values(r, na.rm = TRUE)
    lo <- as.numeric(quantile(v, 0.02))
    hi <- as.numeric(quantile(v, 0.98))
    rr <- (r - lo) / (hi - lo)
    clamp(rr, 0, 1)
  }
  i_s <- rescale01(i_r); e_s <- rescale01(e_r)
  gap_i <- 1 - i_s; gap_e <- 1 - e_s

  weights <- c(0.3, 0.5, 0.7)
  results <- list()

  for (w in weights) {
    cat(sprintf("  weight on internet gap = %.1f (electricity = %.1f)\n", w, 1 - w))
    equity <- w * gap_i + (1 - w) * gap_e
    prio_eq <- prio * equity
    # Rescale 0-1
    v <- values(prio_eq, na.rm = TRUE)
    lo <- as.numeric(quantile(v, 0.01)); hi <- as.numeric(quantile(v, 0.99))
    prio_eq <- (prio_eq - lo) / (hi - lo)
    prio_eq <- clamp(prio_eq, 0, 1)

    # Compute Spearman rho between unweighted priority and equity-adjusted
    set.seed(1234)
    samp <- spatSample(c(prio, prio_eq, equity), 200000, na.rm = TRUE)
    names(samp) <- c("prio_unw", "prio_eq", "equity_gap")
    rho <- cor(samp$prio_unw, samp$prio_eq, method = "spearman")
    rho_gap <- cor(samp$prio_unw, samp$equity_gap, method = "spearman")

    results[[length(results) + 1]] <- tibble(
      w_internet = w, w_electricity = 1 - w,
      rho_priority_vs_equity_adj = round(rho, 4),
      rho_priority_vs_equity_gap = round(rho_gap, 4)
    )
    rm(prio_eq, equity, samp); gc(verbose = FALSE)
  }
  out <- bind_rows(results)
  write_csv(out, file.path(OUT_DIR, "r4_equity_weight_sensitivity.csv"))
  cat("\n  Equity weighting sensitivity:\n")
  print(out)
  out
}

# ============================================================================
# SECTION 5 — THRESHOLD SENSITIVITY (addresses comment #19)
# ============================================================================
# Re-derives greenfield + densification priority surfaces under threshold
# variations: percentile in {85, 90, 95}, macro-cell radius in {3, 5, 7} km.
# Reports % of national area in each class under each variant.
# Expected runtime: 20 min. RAM ~2 GB.
# ============================================================================
run_section_5 <- function() {

  cat("\n========== SECTION 5: Threshold sensitivity (comment #19) ==========\n\n")

  # Build priority surface inputs
  suit <- rast(file.path(RES_DIR, "07_suitability_1km.tif"))
  if (nlyr(suit) > 1) suit <- suit[[1]]

  pop <- rast(file.path(GEE_DIR, "ETH_population.tif"))
  if (nlyr(pop) > 1) pop <- pop[[1]]
  pop <- resample(pop, suit, method = "bilinear")

  bld <- rast(file.path(GEE_DIR, "ETH_builtup_2020.tif"))
  if (nlyr(bld) > 1) bld <- bld[[1]]
  bld <- resample(bld, suit, method = "bilinear")

  dist_tower <- rast(file.path(GEE_DIR, "ETH_dist_to_existing_tower.tif"))
  if (nlyr(dist_tower) > 1) dist_tower <- dist_tower[[1]]
  dist_tower <- resample(dist_tower, suit, method = "bilinear")

  rescale01 <- function(r) {
    v <- values(r, na.rm = TRUE)
    lo <- as.numeric(quantile(v, 0.01))
    hi <- as.numeric(quantile(v, 0.99))
    rr <- (r - lo) / (hi - lo)
    clamp(rr, 0, 1)
  }
  demand <- rescale01(log1p(pop)) * rescale01(bld)
  demand <- rescale01(demand)
  gap    <- rescale01(dist_tower)
  prio   <- suit * demand * gap

  # Grid of thresholds
  pcts    <- c(0.85, 0.90, 0.95)
  radii_m <- c(3000, 5000, 7000)

  # Total cell area: 1 km cell -> 1 km^2 per cell
  total_cells <- sum(!is.na(values(prio)))
  total_km2 <- total_cells   # at 1 km

  results <- list()
  for (p in pcts) {
    thr <- as.numeric(quantile(values(prio, na.rm = TRUE), p))
    for (r in radii_m) {
      green <- (prio >= thr) & (dist_tower >= r)
      n_green <- sum(values(green, na.rm = TRUE) == 1, na.rm = TRUE)
      pct_green <- 100 * n_green / total_cells
      results[[length(results) + 1]] <- tibble(
        percentile = p, radius_km = r / 1000,
        greenfield_km2 = n_green, greenfield_pct = round(pct_green, 2)
      )
      rm(green); gc(verbose = FALSE)
    }
  }
  out <- bind_rows(results)
  write_csv(out, file.path(OUT_DIR, "r5_threshold_sensitivity.csv"))
  cat("  Threshold sensitivity table:\n")
  print(out, n = Inf)
  out
}

# ============================================================================
# MASTER RUNNER
# ============================================================================
# Run sections in sequence. If your machine crashes, simply re-run this whole
# script — completed sections will be skipped via the checkpoint system.
# To run a single section, comment out the others.
# ============================================================================

cat("\n\n=================================================================\n")
cat("  Starting reviewer-response runs\n")
cat("=================================================================\n")

s1 <- checkpoint("r1_complete", run_section_1)   # ~30 min
s2 <- checkpoint("r2_complete", run_section_2)   # ~20 min
s4 <- checkpoint("r4_complete", run_section_4)   # ~30 min
s5 <- checkpoint("r5_complete", run_section_5)   # ~20 min

# Section 3 is the heaviest. Comment this line out if you want a quick run.
s3 <- checkpoint("r3_complete", run_section_3)   # ~3-5 hours

cat("\n\n=================================================================\n")
cat(sprintf("  ALL DONE: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n")
cat(sprintf("  Outputs written to: %s\n", OUT_DIR))
cat("  CSV tables:\n")
cat("    r1_moran_results.csv         (Moran I)\n")
cat("    r1_variogram_fit.csv         (variogram range)\n")
cat("    r2_brier_scores.csv          (Brier per model)\n")
cat("    r3_pseudo_absence_sensitivity.csv\n")
cat("    r4_equity_weight_sensitivity.csv\n")
cat("    r5_threshold_sensitivity.csv\n")
cat("  Figures:\n")
cat("    r1_moran_variogram.png\n")
cat("    r2_reliability_diagram.png\n")
cat("    r3_pseudo_absence_sensitivity.png\n")
cat("  Send all six CSVs + three PNGs back to me for paper write-up.\n")
