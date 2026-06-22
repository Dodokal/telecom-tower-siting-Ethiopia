# ============================================================================
# MINI PATCH: Auto-detect priority surface for Section 4
# ----------------------------------------------------------------------------
# Your file is named "08_priority_score.tif" (no "_1km" suffix).
# This script:
#   1. Auto-detects either filename
#   2. Uses the existing priority surface if found (matches your paper exactly)
#   3. Falls back to rebuilding from suitability if neither exists
#   4. Re-runs only Section 4
#
# To use: save as 10b_section4_fix.R, source it, done in ~10 minutes.
# ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(fs)
})
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.5)

ROOT     <- "K:/ETH TOWERS"
RES_DIR  <- file.path(ROOT, "results")
GEE_DIR  <- file.path(ROOT, "ETH_towers")
OUT_DIR  <- file.path(RES_DIR, "reviewer_response")
CKPT_DIR <- file.path(OUT_DIR, "checkpoints")
dir_create(c(OUT_DIR, CKPT_DIR))

# Delete the old broken checkpoint if it exists
ckpt <- file.path(CKPT_DIR, "r4_complete.rds")
if (file.exists(ckpt)) {
  file.remove(ckpt)
  cat("  Deleted old r4 checkpoint\n")
}

# === AUTO-DETECT PRIORITY SURFACE ===
find_priority <- function() {
  candidates <- c(
    file.path(RES_DIR, "08_priority_score.tif"),       # your file
    file.path(RES_DIR, "08_priority_score_1km.tif"),
    file.path(RES_DIR, "priority_score.tif"),
    file.path(RES_DIR, "08_priority_1km.tif")
  )
  for (p in candidates) {
    if (file.exists(p)) {
      cat(sprintf("  Found priority surface at: %s\n", p))
      return(p)
    }
  }
  cat("  No priority surface found on disk; will rebuild from suitability\n")
  return(NULL)
}

# === LOAD HELPERS ===
load_aligned <- function(path, ref) {
  r <- rast(path); if (nlyr(r) > 1) r <- r[[1]]
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

# ============================================================================
# RUN SECTION 4 — using existing priority surface where possible
# ============================================================================
cat("\n========== SECTION 4 (FINAL): Equity weighting sensitivity ==========\n")
cat(sprintf("  Started: %s\n\n", format(Sys.time(), "%H:%M:%S")))

# Try to use the existing priority surface first
prio_path <- find_priority()

if (!is.null(prio_path)) {
  # Use existing priority surface (preferred — matches your paper exactly)
  prio <- rast(prio_path)
  if (nlyr(prio) > 1) prio <- prio[[1]]
  cat(sprintf("  Loaded priority surface: %d x %d cells\n",
              ncol(prio), nrow(prio)))
  cat("  Range: ", paste(round(range(values(prio, na.rm = TRUE)), 3),
                         collapse = " to "), "\n")

  # Use prio as the reference grid for DHS layers
  ref_grid <- prio
} else {
  # Rebuild from inputs (fallback)
  cat("  Rebuilding priority surface from suitability + demand + gap\n")
  suit <- rast(file.path(RES_DIR, "07_suitability_1km.tif"))
  if (nlyr(suit) > 1) suit <- suit[[1]]
  pop <- load_aligned(file.path(GEE_DIR, "ETH_population.tif"), suit)
  bld <- load_aligned(file.path(GEE_DIR, "ETH_builtup_2020.tif"), suit)
  dist_tower <- load_aligned(file.path(GEE_DIR, "ETH_dist_to_existing_tower.tif"),
                             suit)
  demand <- rescale01(log1p(pop)) * rescale01(bld)
  demand <- rescale01(demand)
  gap    <- rescale01(dist_tower)
  prio   <- suit * demand * gap
  ref_grid <- prio
  rm(suit, pop, bld, dist_tower, demand, gap); gc(verbose = FALSE)
}

# Load DHS layers and align to priority grid
cat("  Loading DHS internet / electricity layers\n")
i_r <- load_aligned(file.path(GEE_DIR, "ETH_dhs_internet_users.tif"), ref_grid)
e_r <- load_aligned(file.path(GEE_DIR, "ETH_dhs_electricity_access.tif"), ref_grid)
i_s <- rescale01(i_r, 0.02, 0.98)
e_s <- rescale01(e_r, 0.02, 0.98)
gap_i <- 1 - i_s
gap_e <- 1 - e_s

weights <- c(0.3, 0.5, 0.7)
results <- list()

for (w in weights) {
  cat(sprintf("\n  weight on internet gap = %.1f (electricity = %.1f)\n",
              w, 1 - w))
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
    w_internet = w,
    w_electricity = 1 - w,
    rho_priority_vs_equity_adj = round(rho_adj, 4),
    rho_priority_vs_equity_gap = round(rho_gap, 4)
  )
  rm(prio_eq, equity, samp); gc(verbose = FALSE)
}

out <- bind_rows(results)
write_csv(out, file.path(OUT_DIR, "r4_equity_weight_sensitivity.csv"))
saveRDS(out, file.path(CKPT_DIR, "r4_complete.rds"))

cat("\n========== EQUITY WEIGHTING SENSITIVITY RESULTS ==========\n")
print(out)

cat("\n\n=================================================================\n")
cat(sprintf("  S4 DONE: %s\n", format(Sys.time(), "%H:%M:%S")))
cat("=================================================================\n")
cat("  Output saved:\n")
cat(sprintf("    %s\n", file.path(OUT_DIR, "r4_equity_weight_sensitivity.csv")))
cat("\n  Send this CSV back to me along with the rest.\n")
