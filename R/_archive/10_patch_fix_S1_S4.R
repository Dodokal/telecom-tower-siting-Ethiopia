# ============================================================================
# PATCH: Fix S1 and S4 errors
# ----------------------------------------------------------------------------
# Problems found:
#   1. training_table.csv is at K:/ETH TOWERS/, not K:/ETH TOWERS/results/
#   2. 08_priority_score_1km.tif was never saved to disk; rebuild on the fly
#
# This patch:
#   - Auto-detects training table location
#   - Rebuilds the priority surface in memory for S4
#   - Forces re-run of S1 and S4 by deleting their checkpoints
#   - Leaves S2, S3, S5 alone (you already completed S2 and S5)
# ============================================================================

# === STEP 0: ESCAPE BROWSER MODE FIRST ===
# If you see "Browse[1]>" at the R prompt, type Q and press Enter before
# running this script. The script will not work inside the debugger.

# === STEP 1: ENVIRONMENT ===
suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(fs)
  library(ggplot2); library(patchwork)
  library(spdep); library(gstat); library(automap)
  library(scales); library(viridis)
})
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.5)

ROOT     <- "K:/ETH TOWERS"
RES_DIR  <- file.path(ROOT, "results")
GEE_DIR  <- file.path(ROOT, "ETH_towers")
OUT_DIR  <- file.path(RES_DIR, "reviewer_response")
FIG_DIR  <- file.path(OUT_DIR, "figures")
CKPT_DIR <- file.path(OUT_DIR, "checkpoints")
dir_create(c(OUT_DIR, FIG_DIR, CKPT_DIR))

theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey30"),
          legend.position = "bottom")
}
save_fig <- function(p, file, w, h)
  ggsave(file, p, width = w, height = h, dpi = 300, bg = "white")

# === STEP 2: AUTO-DETECT TRAINING TABLE LOCATION ===
find_training_table <- function() {
  candidates <- c(
    file.path(ROOT, "training_table.csv"),
    file.path(RES_DIR, "training_table.csv"),
    file.path(ROOT, "results", "training_table.csv")
  )
  for (p in candidates) {
    if (file.exists(p)) {
      cat(sprintf("  Found training table at: %s\n", p))
      return(p)
    }
  }
  stop("training_table.csv not found at any expected location.\n",
       "  Checked:\n    ", paste(candidates, collapse = "\n    "))
}

find_presences <- function() {
  candidates <- c(
    file.path(ROOT, "presence_absence.gpkg"),
    file.path(RES_DIR, "presence_absence.gpkg"),
    file.path(ROOT, "presence_absence.csv"),
    file.path(RES_DIR, "presence_absence.csv")
  )
  for (p in candidates) {
    if (file.exists(p)) {
      cat(sprintf("  Found presence/absence at: %s\n", p))
      return(p)
    }
  }
  stop("presence_absence file not found.")
}

# === STEP 3: DELETE OLD BROKEN CHECKPOINTS ===
broken_ckpts <- c("r1_complete.rds", "r4_complete.rds")
for (f in broken_ckpts) {
  full <- file.path(CKPT_DIR, f)
  if (file.exists(full)) {
    file.remove(full)
    cat(sprintf("  Deleted broken checkpoint: %s\n", f))
  }
}

# ============================================================================
# FIXED SECTION 1 — MORAN I + VARIOGRAM
# ============================================================================
run_section_1_fixed <- function() {

  cat("\n========== FIXED SECTION 1: Moran I + variogram ==========\n\n")

  tt_path <- find_training_table()
  pa_path <- find_presences()

  tt <- read_csv(tt_path, show_col_types = FALSE)
  cat(sprintf("  training_table: %d rows x %d cols\n", nrow(tt), ncol(tt)))

  suit_r <- rast(file.path(RES_DIR, "07_suitability_1km.tif"))
  if (nlyr(suit_r) > 1) suit_r <- suit_r[[1]]

  # Load presences (handle both gpkg and csv)
  if (grepl("\\.gpkg$", pa_path)) {
    pa <- st_read(pa_path, quiet = TRUE)
  } else {
    pa <- read_csv(pa_path, show_col_types = FALSE)
    # Detect coordinate columns
    xcol <- intersect(c("x", "X", "lon", "longitude", "Longitude"), names(pa))[1]
    ycol <- intersect(c("y", "Y", "lat", "latitude", "Latitude"), names(pa))[1]
    if (is.na(xcol) || is.na(ycol))
      stop("Cannot find coordinate columns in presence/absence CSV")
    pa <- st_as_sf(pa, coords = c(xcol, ycol), crs = 4326)
  }

  # Identify class column
  class_col <- intersect(c("class", "Class", "label", "presence", "pres"),
                        names(pa))[1]
  if (is.na(class_col)) {
    if ("pres" %in% names(pa)) pa$class <- ifelse(pa$pres == 1, "presence", "absence")
    else stop("Cannot identify class column in presence/absence data")
  } else {
    if (class_col != "class") pa$class <- pa[[class_col]]
  }

  # Get predictions at points
  pa_pred <- terra::extract(suit_r, vect(st_transform(pa, crs(suit_r))))
  pa$pred <- pa_pred[, 2]
  pa$obs  <- as.integer(grepl("pres|^1$|TRUE", as.character(pa$class),
                              ignore.case = TRUE))
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

  # Variogram
  cat("  Fitting empirical variogram\n")
  pa_sp <- as_Spatial(pa_utm[keep, ])
  pa_sp$resid <- resid
  vg <- variogram(resid ~ 1, data = pa_sp, cutoff = 200000, width = 10000)
  vg_fit <- tryCatch(
    autofitVariogram(resid ~ 1, input_data = pa_sp,
                     model = c("Sph", "Exp", "Gau"))$var_model,
    error = function(e) {
      cat("  [warn] autoKrige failed, manual exponential fit\n")
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

  p_all <- p_moran | p_vg
  save_fig(p_all, file.path(FIG_DIR, "r1_moran_variogram.png"), 11, 4.5)
  write_csv(moran_df, file.path(OUT_DIR, "r1_moran_results.csv"))
  write_csv(tibble(model = vg_fit$model[2],
                   range_km = round(range_m / 1000, 2),
                   psill = vg_fit$psill[2],
                   nugget = vg_fit$psill[1]),
            file.path(OUT_DIR, "r1_variogram_fit.csv"))

  cat("  S1 done. Files saved:\n")
  cat("    r1_moran_results.csv\n    r1_variogram_fit.csv\n")
  cat("    figures/r1_moran_variogram.png\n")

  list(moran = moran_df, vg_range_km = range_m / 1000)
}

# ============================================================================
# FIXED SECTION 4 — EQUITY WEIGHTING SENSITIVITY
# ============================================================================
# Rebuilds the priority surface from inputs you already have.
# ============================================================================
run_section_4_fixed <- function() {

  cat("\n========== FIXED SECTION 4: Equity weighting sensitivity ==========\n\n")

  # Load the LightGBM suitability raster (this DOES exist on your disk)
  suit_path <- file.path(RES_DIR, "07_suitability_1km.tif")
  if (!file.exists(suit_path))
    stop("07_suitability_1km.tif not found at ", suit_path)

  suit <- rast(suit_path)
  if (nlyr(suit) > 1) suit <- suit[[1]]
  cat(sprintf("  Loaded suitability surface: %d x %d cells\n",
              ncol(suit), nrow(suit)))

  # Helper: load single-layer raster, resample to suit grid
  load_aligned <- function(path) {
    r <- rast(path); if (nlyr(r) > 1) r <- r[[1]]
    if (!compareGeom(r, suit, stopOnError = FALSE)) {
      r <- resample(r, suit, method = "bilinear", threads = TRUE)
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

  # REBUILD the priority surface in memory
  cat("  Rebuilding priority surface from inputs\n")
  pop <- load_aligned(file.path(GEE_DIR, "ETH_population.tif"))
  bld <- load_aligned(file.path(GEE_DIR, "ETH_builtup_2020.tif"))
  dist_tower <- load_aligned(file.path(GEE_DIR, "ETH_dist_to_existing_tower.tif"))

  demand <- rescale01(log1p(pop)) * rescale01(bld)
  demand <- rescale01(demand)
  gap    <- rescale01(dist_tower)
  prio   <- suit * demand * gap

  cat("  Priority surface rebuilt\n")
  rm(pop, bld); gc(verbose = FALSE)

  # Load DHS layers
  i_r <- load_aligned(file.path(GEE_DIR, "ETH_dhs_internet_users.tif"))
  e_r <- load_aligned(file.path(GEE_DIR, "ETH_dhs_electricity_access.tif"))
  i_s <- rescale01(i_r, 0.02, 0.98); e_s <- rescale01(e_r, 0.02, 0.98)
  gap_i <- 1 - i_s; gap_e <- 1 - e_s

  weights <- c(0.3, 0.5, 0.7)
  results <- list()

  for (w in weights) {
    cat(sprintf("  weight on internet gap = %.1f (electricity = %.1f)\n",
                w, 1 - w))
    equity <- w * gap_i + (1 - w) * gap_e
    prio_eq <- prio * equity
    prio_eq <- rescale01(prio_eq)

    set.seed(1234)
    samp <- spatSample(c(prio, prio_eq, equity), 200000,
                       na.rm = TRUE, as.df = TRUE)
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
  cat("\n  S4 done. File saved:\n")
  cat("    r4_equity_weight_sensitivity.csv\n")
  out
}

# ============================================================================
# RUN BOTH FIXES
# ============================================================================
cat("\n=================================================================\n")
cat("  Running fixed S1 and S4\n")
cat(sprintf("  Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n")

# Save checkpoints so we can resume / skip
s1_out <- run_section_1_fixed()
saveRDS(s1_out, file.path(CKPT_DIR, "r1_complete.rds"))

s4_out <- run_section_4_fixed()
saveRDS(s4_out, file.path(CKPT_DIR, "r4_complete.rds"))

cat("\n\n=================================================================\n")
cat(sprintf("  PATCH DONE: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=================================================================\n")
cat("  All four light sections (S1, S2, S4, S5) are now complete.\n")
cat("  Send me back:\n")
cat("    r1_moran_results.csv\n")
cat("    r1_variogram_fit.csv\n")
cat("    r2_brier_scores.csv\n")
cat("    r4_equity_weight_sensitivity.csv\n")
cat("    r5_threshold_sensitivity.csv\n")
cat("    figures/r1_moran_variogram.png\n")
cat("    figures/r2_reliability_diagram.png\n")
