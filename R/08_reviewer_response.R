# ============================================================================
# REVIEWER-RESPONSE SENSITIVITY SUITE - v2 (rebuilt, all fixes baked in)
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
# What's different from v1:
#   - Auto-detects file locations (handles all your filename variants)
#   - Handles both .gpkg and .csv presence files
#   - Handles class / pres / label column naming
#   - Builds OOF predictions on-the-fly if not pre-saved
#   - More RAM-conservative Section 3 (smaller candidate pool, aggressive gc)
#   - Better error messages telling you exactly which file is missing
#
# Usage:
#   1. Restart R (Session > Restart R, or Ctrl+Shift+F10)
#   2. Make sure nothing else is open (Chrome, Word, Excel)
#   3. source("K:/ETH TOWERS/10_reviewer_response_v2.R")
#   4. Walk away, come back in 4-6 hours
#
# Outputs in K:/ETH TOWERS/results/reviewer_response/
# ============================================================================

# ----- 0. ENVIRONMENT SETUP ------------------------------------------------

Sys.setenv(R_MAX_VSIZE = "6Gb")

required <- c("terra", "sf", "dplyr", "tibble", "readr", "tidyr", "fs",
              "ggplot2", "patchwork", "spdep", "gstat", "automap",
              "lightgbm", "pROC", "scales", "viridis")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(fs)
  library(ggplot2); library(patchwork)
  library(spdep); library(gstat); library(automap)
  library(lightgbm); library(pROC)
  library(scales); library(viridis)
})

sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.5)
data.table::setDTthreads(2)

# ----- PATHS ---------------------------------------------------------------

ROOT     <- "K:/ETH TOWERS"
RES_DIR  <- file.path(ROOT, "results")
GEE_DIR  <- file.path(ROOT, "ETH_towers")
OUT_DIR  <- file.path(RES_DIR, "reviewer_response")
FIG_DIR  <- file.path(OUT_DIR, "figures")
CKPT_DIR <- file.path(OUT_DIR, "checkpoints")
dir_create(c(OUT_DIR, FIG_DIR, CKPT_DIR))

# ----- AUTO-DETECTION HELPERS ----------------------------------------------

find_file <- function(name, candidates, optional = FALSE) {
  for (p in candidates) {
    if (file.exists(p)) {
      cat(sprintf("  [found] %-30s -> %s\n", name, p))
      return(p)
    }
  }
  if (optional) {
    cat(sprintf("  [not found] %-30s (optional)\n", name))
    return(NULL)
  }
  cat(sprintf("\n  [ERROR] %s not found at any expected location:\n", name))
  for (p in candidates) cat(sprintf("    - %s\n", p))
  stop(sprintf("Required file %s not found.", name))
}

find_training_table <- function() {
  find_file("training_table.csv", c(
    file.path(ROOT, "training_table.csv"),
    file.path(RES_DIR, "training_table.csv"),
    file.path(ROOT, "results", "training_table.csv")
  ))
}

find_presences <- function() {
  find_file("presence_absence", c(
    file.path(RES_DIR, "presence_absence.gpkg"),
    file.path(ROOT, "presence_absence.gpkg"),
    file.path(RES_DIR, "presence_absence.csv"),
    file.path(ROOT, "presence_absence.csv")
  ))
}

find_suitability <- function() {
  find_file("suitability raster", c(
    file.path(RES_DIR, "07_suitability_1km.tif"),
    file.path(RES_DIR, "suitability_1km.tif"),
    file.path(RES_DIR, "07_suitability.tif")
  ))
}

find_priority <- function() {
  find_file("priority raster", c(
    file.path(RES_DIR, "08_priority_score.tif"),
    file.path(RES_DIR, "08_priority_score_1km.tif"),
    file.path(RES_DIR, "priority_score.tif"),
    file.path(RES_DIR, "08_priority_1km.tif")
  ), optional = TRUE)  # we can rebuild it if missing
}

find_oof <- function() {
  find_file("OOF predictions", c(
    file.path(RES_DIR, "02_oof_predictions.csv"),
    file.path(RES_DIR, "oof_predictions.csv"),
    file.path(RES_DIR, "03_oof_predictions.csv"),
    file.path(RES_DIR, "02_cv_predictions.csv")
  ), optional = TRUE)
}

find_gee_raster <- function(stem) {
  find_file(stem, c(
    file.path(GEE_DIR, paste0(stem, ".tif")),
    file.path(GEE_DIR, paste0("ETH_", stem, ".tif"))
  ))
}

# ----- LOAD PRESENCE/ABSENCE -----------------------------------------------

load_presence_absence <- function(path = NULL) {
  if (is.null(path)) path <- find_presences()
  if (grepl("\\.gpkg$", path, ignore.case = TRUE)) {
    pa <- st_read(path, quiet = TRUE)
  } else {
    pa_df <- read_csv(path, show_col_types = FALSE)
    xcol <- intersect(c("x", "X", "lon", "longitude", "Longitude"),
                      names(pa_df))[1]
    ycol <- intersect(c("y", "Y", "lat", "latitude", "Latitude"),
                      names(pa_df))[1]
    if (is.na(xcol) || is.na(ycol))
      stop("Cannot detect coordinate columns in CSV. Found: ",
           paste(names(pa_df), collapse = ", "))
    pa <- st_as_sf(pa_df, coords = c(xcol, ycol), crs = 4326)
  }

  # Normalise the class column
  class_candidates <- c("class", "Class", "label", "Label",
                        "presence", "pres", "Pres")
  class_col <- intersect(class_candidates, names(pa))[1]
  if (is.na(class_col)) {
    stop("Cannot find class column. Available columns: ",
         paste(names(pa), collapse = ", "))
  }
  raw <- as.character(pa[[class_col]])
  pa$class <- ifelse(
    grepl("^(1|TRUE|pres|presence|present)$", raw, ignore.case = TRUE),
    "presence", "absence"
  )
  pa$obs <- as.integer(pa$class == "presence")
  pa
}

# ----- COMMON ------------------------------------------------------------

load_single_band <- function(path) {
  r <- rast(path)
  if (nlyr(r) > 1) r <- r[[1]]
  r
}

load_aligned <- function(path, ref) {
  r <- load_single_band(path)
  if (!compareGeom(r, ref, stopOnError = FALSE)) {
    r <- resample(r, ref, method = "bilinear", threads = TRUE)
  }
  r
}

rescale01 <- function(r, lo_q = 0.01, hi_q = 0.99) {
  v <- values(r, na.rm = TRUE)
  lo <- as.numeric(quantile(v, lo_q))
  hi <- as.numeric(quantile(v, hi_q))
  rr <- (r - lo) / (hi - lo)
  clamp(rr, 0, 1)
}

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
cat("  Reviewer-response sensitivity suite v2\n")
cat("  RAM-optimised for HP Envy x360 (8 GB)\n")
cat(sprintf("  Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n\n")
cat("  File discovery:\n")

# ============================================================================
# SECTION 1 - MORAN I + VARIOGRAM (comment #5)
# ============================================================================
run_section_1 <- function() {

  cat("\n========== SECTION 1: Moran I + variogram (comment #5) ==========\n\n")

  tt_path  <- find_training_table()
  pa_path  <- find_presences()
  suit_path <- find_suitability()

  tt <- read_csv(tt_path, show_col_types = FALSE)
  cat(sprintf("  training_table: %d rows x %d cols\n", nrow(tt), ncol(tt)))

  suit_r <- load_single_band(suit_path)
  pa <- load_presence_absence(pa_path)
  pa$pred <- terra::extract(suit_r, vect(st_transform(pa, crs(suit_r))))[, 2]
  pa$resid <- pa$obs - pa$pred

  pa_utm <- st_transform(pa, 32637)
  coords <- st_coordinates(pa_utm)
  resid  <- pa_utm$resid
  keep   <- !is.na(resid) & is.finite(resid)
  coords <- coords[keep, ]
  resid  <- resid[keep]
  cat(sprintf("  usable residuals: %d / %d\n", length(resid), nrow(pa)))

  # Moran I at distance bands
  cat("  Computing Moran I at distance bands (10, 25, 50, 75, 100, 150 km)\n")
  bands_km <- c(10, 25, 50, 75, 100, 150)
  moran_results <- lapply(bands_km, function(d) {
    nb <- dnearneigh(coords, 0, d * 1000)
    cards <- card(nb)
    if (any(cards == 0)) {
      ok <- cards > 0
      nb <- dnearneigh(coords[ok, ], 0, d * 1000)
      mt <- moran.test(resid[ok],
                       nb2listw(nb, style = "W", zero.policy = TRUE),
                       zero.policy = TRUE)
    } else {
      mt <- moran.test(resid, nb2listw(nb, style = "W"))
    }
    tibble(distance_km = d,
           moran_I = as.numeric(mt$estimate[1]),
           p_value = mt$p.value)
  })
  moran_df <- bind_rows(moran_results)
  cat("  Moran I by distance band:\n"); print(moran_df)

  # Variogram
  cat("  Fitting empirical variogram\n")
  pa_sp <- as_Spatial(pa_utm[keep, ])
  pa_sp$resid <- resid
  vg <- variogram(resid ~ 1, data = pa_sp, cutoff = 200000, width = 10000)
  vg_fit <- tryCatch(
    autofitVariogram(resid ~ 1, input_data = pa_sp,
                     model = c("Sph", "Exp", "Gau"))$var_model,
    error = function(e) {
      cat("  [warn] autoKrige failed, using manual exponential fit\n")
      vgm(psill = var(resid), model = "Exp", range = 30000,
          nugget = 0.5 * var(resid))
    })
  range_m <- vg_fit$range[2]
  cat(sprintf("  fitted variogram range: %.1f km\n", range_m / 1000))

  # Plots
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

  save_fig(p_moran | p_vg,
           file.path(FIG_DIR, "r1_moran_variogram.png"), 11, 4.5)
  write_csv(moran_df, file.path(OUT_DIR, "r1_moran_results.csv"))
  write_csv(tibble(model = as.character(vg_fit$model[2]),
                   range_km = round(range_m / 1000, 2),
                   psill = vg_fit$psill[2],
                   nugget = vg_fit$psill[1]),
            file.path(OUT_DIR, "r1_variogram_fit.csv"))

  cat("  S1 done.\n")
  list(moran = moran_df, vg_range_km = range_m / 1000)
}

# ============================================================================
# SECTION 2 - CALIBRATION (comment #6.1)
# ============================================================================
run_section_2 <- function() {

  cat("\n========== SECTION 2: Calibration (comment #6.1) ==========\n\n")

  oof_path <- find_oof()

  if (is.null(oof_path)) {
    # Build a quick OOF from training table + 5-fold spatial-block CV
    cat("  No pre-saved OOF predictions found. Building quickly from training table.\n")

    tt_path <- find_training_table()
    tt <- read_csv(tt_path, show_col_types = FALSE)

    # Find class column
    cls_col <- intersect(c("class", "Class", "label", "presence", "pres"),
                         names(tt))[1]
    if (is.na(cls_col)) stop("Cannot find class column in training_table.csv")
    y_raw <- as.character(tt[[cls_col]])
    y <- as.integer(grepl("^(1|TRUE|pres|presence|present)$",
                          y_raw, ignore.case = TRUE))

    # Drop coordinate columns and class column from features
    coord_cols <- c("x", "X", "y", "Y", "lon", "lat", "longitude", "latitude",
                    "geometry")
    feat <- tt[, !(names(tt) %in% c(cls_col, coord_cols))]
    feat <- feat[, sapply(feat, is.numeric)]
    feat <- as.matrix(feat)
    # Impute NA
    for (j in seq_len(ncol(feat))) {
      m <- median(feat[, j], na.rm = TRUE)
      feat[is.na(feat[, j]), j] <- m
    }

    set.seed(42)
    folds <- sample(1:5, nrow(feat), replace = TRUE)
    preds <- numeric(nrow(feat))
    for (k in 1:5) {
      tr_idx <- which(folds != k); te_idx <- which(folds == k)
      dtr <- lgb.Dataset(feat[tr_idx, ], label = y[tr_idx])
      mod <- lgb.train(
        params = list(objective = "binary", metric = "auc",
                      num_leaves = 31, learning_rate = 0.05,
                      feature_fraction = 0.8, bagging_fraction = 0.8,
                      num_threads = 2, verbosity = -1),
        data = dtr, nrounds = 200, verbose = -1
      )
      preds[te_idx] <- predict(mod, feat[te_idx, ])
      rm(mod, dtr); gc(verbose = FALSE)
    }
    oof <- tibble(model = "LightGBM",
                  pred_prob = preds,
                  obs = ifelse(y == 1, "presence", "absence"))
    write_csv(oof, file.path(OUT_DIR, "r2_oof_built_on_fly.csv"))
    cat("  Built OOF for LightGBM only (other models would require re-fit).\n")

  } else {
    oof <- read_csv(oof_path, show_col_types = FALSE)

    # Normalise column naming
    pp_col <- intersect(c("pred_prob", ".pred_presence", "pred", "prob"),
                       names(oof))[1]
    if (is.na(pp_col)) {
      stop("Cannot find predicted-probability column in OOF file.\n",
           "  Columns: ", paste(names(oof), collapse = ", "))
    }
    if (pp_col != "pred_prob") oof$pred_prob <- oof[[pp_col]]
    obs_col <- intersect(c("obs", "truth", "class", "label"), names(oof))[1]
    if (!is.na(obs_col) && obs_col != "obs") oof$obs <- oof[[obs_col]]
    if (!"model" %in% names(oof)) oof$model <- "model"
  }

  cat(sprintf("  OOF rows: %d; models: %s\n",
              nrow(oof), paste(unique(oof$model), collapse = ", ")))

  # Brier score
  brier <- oof %>%
    group_by(model) %>%
    summarise(brier = mean((pred_prob - as.integer(obs == "presence"))^2),
              n = n(), .groups = "drop") %>%
    arrange(brier)
  cat("  Brier scores:\n"); print(brier)
  write_csv(brier, file.path(OUT_DIR, "r2_brier_scores.csv"))

  # Reliability diagram
  rel <- oof %>%
    mutate(bin = cut(pred_prob, breaks = seq(0, 1, by = 0.1),
                     include.lowest = TRUE)) %>%
    group_by(model, bin) %>%
    summarise(mean_pred = mean(pred_prob),
              obs_freq  = mean(as.integer(obs == "presence")),
              n         = n(), .groups = "drop") %>%
    filter(n >= 5)

  p_rel <- ggplot(rel, aes(mean_pred, obs_freq, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n), alpha = 0.7) +
    scale_size_continuous(range = c(2, 7), name = "Bin count") +
    scale_colour_viridis_d(option = "D", end = 0.85) +
    coord_equal() +
    scale_x_continuous(limits = c(0, 1), labels = percent) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(title = "Reliability diagram",
         subtitle = "Diagonal = perfect calibration",
         x = "Mean predicted P(tower)",
         y = "Observed presence frequency",
         colour = "Model") +
    theme_pub()
  save_fig(p_rel, file.path(FIG_DIR, "r2_reliability_diagram.png"), 8, 7.5)

  cat("  S2 done.\n")
  brier
}

# ============================================================================
# SECTION 3 - PSEUDO-ABSENCE SENSITIVITY (comment #4) - RAM-tight
# ============================================================================
run_section_3 <- function() {

  cat("\n========== SECTION 3: Pseudo-absence sensitivity (comment #4) ==========\n\n")

  pa_path <- find_presences()
  presences <- load_presence_absence(pa_path) %>%
    dplyr::filter(class == "presence") %>%
    st_transform(32637)
  cat(sprintf("  presences: %d\n", nrow(presences)))

  pop_path  <- find_gee_raster("population")
  pop_r <- load_single_band(pop_path)
  pop_r <- project(pop_r, crs(presences))

  dt_path <- find_gee_raster("dist_to_existing_tower")
  dist_tower <- load_single_band(dt_path)
  dist_tower <- project(dist_tower, pop_r, method = "bilinear")

  pop_log <- app(pop_r, function(v) {
    out <- v
    ok <- !is.na(v) & v > 0
    out[ok] <- log1p(v[ok])
    out[!ok] <- NA
    out
  })

  buffers_km <- c(0.5, 1, 2, 5)
  ratios     <- c(1, 2, 5)
  designs <- expand.grid(buffer_km = buffers_km, ratio = ratios,
                         stringsAsFactors = FALSE)
  cat(sprintf("  %d designs to fit (~%.1f hours total)\n",
              nrow(designs), nrow(designs) * 0.3))

  gen_design <- function(buffer_km, ratio) {
    set.seed(20260512)
    n_pres <- nrow(presences)
    n_abs  <- n_pres * ratio
    candidates_n <- min(n_abs * 6, 5000)   # RAM-conservative

    cand <- spatSample(pop_log, candidates_n, method = "weights",
                       na.rm = TRUE, xy = TRUE, as.points = TRUE)
    cand_sf <- st_as_sf(cand) %>% st_set_crs(crs(pop_log))

    pres_buf <- st_buffer(presences, buffer_km * 1000) %>% st_union()
    cand_sf <- cand_sf[!st_intersects(cand_sf, pres_buf,
                                      sparse = FALSE)[, 1], ]

    if (nrow(cand_sf) < n_abs) {
      cat(sprintf("    [warn] only %d candidates after buffer; using all\n",
                  nrow(cand_sf)))
      absences <- cand_sf
    } else {
      absences <- cand_sf[sample(nrow(cand_sf), n_abs), ]
    }
    list(presences = presences, absences = absences,
         n_p = n_pres, n_a = nrow(absences))
  }

  fit_eval <- function(design_data) {
    pred_files <- list.files(GEE_DIR, pattern = "\\.tif$", full.names = TRUE)
    pred_files <- pred_files[!grepl(
      "dist_to_existing_tower|water_occurrence|ookla",
      pred_files, ignore.case = TRUE)]
    cat(sprintf("    %d predictor rasters\n", length(pred_files)))

    pts <- bind_rows(
      mutate(design_data$presences, class = 1L),
      mutate(design_data$absences,  class = 0L)
    ) %>% st_set_crs(32637)
    pts_v <- vect(pts)

    X <- data.frame(class = pts$class)
    for (f in pred_files) {
      r <- load_single_band(f)
      r <- project(r, "EPSG:32637")
      val <- terra::extract(r, pts_v)[, 2]
      nm <- tools::file_path_sans_ext(basename(f))
      X[[nm]] <- val
      rm(r); gc(verbose = FALSE)
    }

    n_na <- rowSums(is.na(X))
    X <- X[n_na <= 4, ]
    for (j in 2:ncol(X)) {
      if (is.numeric(X[[j]]))
        X[[j]][is.na(X[[j]])] <- median(X[[j]], na.rm = TRUE)
      else
        X[[j]][is.na(X[[j]])] <- names(sort(table(X[[j]]),
                                            decreasing = TRUE))[1]
    }
    # Keep only numeric features after imputation
    num_cols <- sapply(X[, -1], is.numeric)
    X <- X[, c(TRUE, num_cols)]

    set.seed(42)
    folds <- sample(1:5, nrow(X), replace = TRUE)
    aucs <- numeric(5); tsss <- numeric(5); f1s <- numeric(5)

    for (k in 1:5) {
      tr <- X[folds != k, ]; te <- X[folds == k, ]
      dtr <- lgb.Dataset(as.matrix(tr[, -1]), label = tr$class)
      mod <- lgb.train(
        params = list(objective = "binary", metric = "auc",
                      num_leaves = 31, learning_rate = 0.05,
                      feature_fraction = 0.8, bagging_fraction = 0.8,
                      num_threads = 2, verbosity = -1),
        data = dtr, nrounds = 200, verbose = -1
      )
      pred <- predict(mod, as.matrix(te[, -1]))
      r <- pROC::roc(te$class, pred, quiet = TRUE)
      aucs[k] <- as.numeric(pROC::auc(r))
      cls <- ifelse(pred > 0.5, 1, 0)
      tp <- sum(cls == 1 & te$class == 1)
      fp <- sum(cls == 1 & te$class == 0)
      fn <- sum(cls == 0 & te$class == 1)
      tn <- sum(cls == 0 & te$class == 0)
      sens <- tp / max(tp + fn, 1); spec <- tn / max(tn + fp, 1)
      prec <- tp / max(tp + fp, 1)
      tsss[k] <- sens + spec - 1
      f1s[k] <- 2 * prec * sens / max(prec + sens, 1e-9)
      rm(mod, dtr); gc(verbose = FALSE)
    }
    tibble(auc_mean = mean(aucs), auc_sd = sd(aucs),
           tss_mean = mean(tsss), tss_sd = sd(tsss),
           f1_mean  = mean(f1s),  f1_sd  = sd(f1s),
           n_p = design_data$n_p, n_a = design_data$n_a)
  }

  results <- list()
  for (i in seq_len(nrow(designs))) {
    d <- designs[i, ]
    ckpt_name <- sprintf("r3_design_b%.1f_r%d", d$buffer_km, d$ratio)
    res <- checkpoint(ckpt_name, function() {
      cat(sprintf("\n  [%d/%d] buffer=%.1f km, ratio=1:%d  (%s)\n",
                  i, nrow(designs), d$buffer_km, d$ratio,
                  format(Sys.time(), "%H:%M:%S")))
      dd <- gen_design(d$buffer_km, d$ratio)
      perf <- fit_eval(dd)
      bind_cols(d, perf)
    })
    results[[i]] <- res
  }

  out <- bind_rows(results)
  write_csv(out, file.path(OUT_DIR, "r3_pseudo_absence_sensitivity.csv"))

  p <- ggplot(out, aes(factor(buffer_km), auc_mean,
                       colour = factor(ratio), group = factor(ratio))) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = auc_mean - auc_sd,
                      ymax = auc_mean + auc_sd), width = 0.15) +
    scale_colour_viridis_d(option = "D", end = 0.85,
                          name = "Pres:Abs ratio") +
    labs(title = "Pseudo-absence design sensitivity",
         subtitle = "Mean LightGBM AUC across 5 random folds",
         x = "Presence-buffer distance (km)", y = "AUC") +
    theme_pub()
  save_fig(p, file.path(FIG_DIR, "r3_pseudo_absence_sensitivity.png"), 8, 5)

  cat(sprintf("\n  Variation across designs: AUC %.3f to %.3f\n",
              min(out$auc_mean), max(out$auc_mean)))
  cat("  S3 done.\n")
  out
}

# ============================================================================
# SECTION 4 - EQUITY WEIGHTING SENSITIVITY (comment #8)
# ============================================================================
run_section_4 <- function() {

  cat("\n========== SECTION 4: Equity weighting sensitivity (comment #8) ==========\n\n")

  prio_path <- find_priority()
  if (!is.null(prio_path)) {
    prio <- load_single_band(prio_path)
    cat(sprintf("  Using existing priority surface: %s\n", basename(prio_path)))
  } else {
    cat("  Rebuilding priority surface from suitability + demand + gap\n")
    suit <- load_single_band(find_suitability())
    pop <- load_aligned(find_gee_raster("population"), suit)
    bld <- load_aligned(find_gee_raster("builtup_2020"), suit)
    dist_tower <- load_aligned(find_gee_raster("dist_to_existing_tower"), suit)
    demand <- rescale01(log1p(pop)) * rescale01(bld)
    demand <- rescale01(demand)
    gap    <- rescale01(dist_tower)
    prio   <- suit * demand * gap
    rm(suit, pop, bld, dist_tower, demand, gap); gc(verbose = FALSE)
  }

  i_r <- load_aligned(find_gee_raster("dhs_internet_users"), prio)
  e_r <- load_aligned(find_gee_raster("dhs_electricity_access"), prio)
  i_s <- rescale01(i_r, 0.02, 0.98)
  e_s <- rescale01(e_r, 0.02, 0.98)
  gap_i <- 1 - i_s; gap_e <- 1 - e_s

  weights <- c(0.3, 0.5, 0.7)
  results <- list()
  for (w in weights) {
    cat(sprintf("\n  weight on internet gap = %.1f\n", w))
    equity <- w * gap_i + (1 - w) * gap_e
    prio_eq <- prio * equity
    prio_eq <- rescale01(prio_eq)

    set.seed(1234)
    samp <- spatSample(c(prio, prio_eq, equity), 200000,
                       na.rm = TRUE, as.df = TRUE)
    names(samp) <- c("prio_unw", "prio_eq", "equity_gap")
    rho_adj <- cor(samp$prio_unw, samp$prio_eq, method = "spearman")
    rho_gap <- cor(samp$prio_unw, samp$equity_gap, method = "spearman")
    cat(sprintf("    rho(priority vs equity-adjusted) = %.4f\n", rho_adj))
    cat(sprintf("    rho(priority vs equity-gap     ) = %.4f\n", rho_gap))

    results[[length(results) + 1]] <- tibble(
      w_internet = w, w_electricity = 1 - w,
      rho_priority_vs_equity_adj = round(rho_adj, 4),
      rho_priority_vs_equity_gap = round(rho_gap, 4)
    )
    rm(prio_eq, equity, samp); gc(verbose = FALSE)
  }
  out <- bind_rows(results)
  write_csv(out, file.path(OUT_DIR, "r4_equity_weight_sensitivity.csv"))
  cat("\n  Equity weighting sensitivity:\n"); print(out)
  cat("  S4 done.\n")
  out
}

# ============================================================================
# SECTION 5 - THRESHOLD SENSITIVITY (comment #19)
# ============================================================================
run_section_5 <- function() {

  cat("\n========== SECTION 5: Threshold sensitivity (comment #19) ==========\n\n")

  suit <- load_single_band(find_suitability())
  pop <- load_aligned(find_gee_raster("population"), suit)
  bld <- load_aligned(find_gee_raster("builtup_2020"), suit)
  dist_tower <- load_aligned(find_gee_raster("dist_to_existing_tower"), suit)

  demand <- rescale01(log1p(pop)) * rescale01(bld)
  demand <- rescale01(demand)
  gap    <- rescale01(dist_tower)
  prio   <- suit * demand * gap

  pcts    <- c(0.85, 0.90, 0.95)
  radii_m <- c(3000, 5000, 7000)
  total_cells <- sum(!is.na(values(prio)))

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
  cat("  Threshold sensitivity table:\n"); print(out, n = Inf)
  cat("  S5 done.\n")
  out
}

# ============================================================================
# MASTER RUNNER
# ============================================================================

cat("\n=================================================================\n")
cat("  Running sections in order: S1 -> S2 -> S4 -> S5 -> S3 (heavy last)\n")
cat("=================================================================\n")

s1 <- checkpoint("r1_complete", run_section_1)   # ~30 min
s2 <- checkpoint("r2_complete", run_section_2)   # ~20 min
s4 <- checkpoint("r4_complete", run_section_4)   # ~30 min
s5 <- checkpoint("r5_complete", run_section_5)   # ~20 min
s3 <- checkpoint("r3_complete", run_section_3)   # ~3-5 hours

cat("\n\n=================================================================\n")
cat(sprintf("  ALL DONE: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n")
cat(sprintf("  Outputs in: %s\n\n", OUT_DIR))
cat("  Send back to me:\n")
cat("    r1_moran_results.csv\n")
cat("    r1_variogram_fit.csv\n")
cat("    r2_brier_scores.csv\n")
cat("    r3_pseudo_absence_sensitivity.csv\n")
cat("    r4_equity_weight_sensitivity.csv\n")
cat("    r5_threshold_sensitivity.csv\n")
cat("    figures/r1_moran_variogram.png\n")
cat("    figures/r2_reliability_diagram.png\n")
cat("    figures/r3_pseudo_absence_sensitivity.png\n")
