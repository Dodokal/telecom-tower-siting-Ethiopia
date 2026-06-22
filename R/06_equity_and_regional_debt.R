# ============================================================================
# Ethiopia Telecom Tower Siting - Equity & Regional Deployment Debt Analysis
# ----------------------------------------------------------------------------
# Adds two distinctive analyses to the manuscript:
#
#   (1) Digital-equity-adjusted priority surface
#       Combines the modelled priority score with two equity weights:
#         * gendered/general internet access gap from DHS 2016
#           (1 - normalised mode internet/mobile users)
#         * electricity access gap (1 - normalised electricity access)
#       The rationale: pixels with high siting suitability AND strong
#       under-access on these dimensions are the highest-impact
#       deployments for closing the digital divide.
#
#   (2) Regional-state deployment debt
#       Aggregates greenfield + densification priority counts per ADM1
#       region, normalises by area and by WorldPop population, and
#       ranks Ethiopia's regional states by their unmet-deployment
#       intensity. Outputs a tabular summary plus a ranked bar chart.
#
# Inputs:
#   K:/ETH TOWERS/results/07_suitability_1km.tif
#   K:/ETH TOWERS/results/08_priority_score.tif
#   K:/ETH TOWERS/results/08_priority_greenfield.tif
#   K:/ETH TOWERS/results/08_priority_densification.tif
#   K:/ETH TOWERS/ETH_towers/ETH_dhs_internet_users.tif
#   K:/ETH TOWERS/ETH_towers/ETH_dhs_electricity_access.tif
#   K:/ETH TOWERS/ETH_towers/ETH_population.tif
#   K:/ETH TOWERS/eth_admin_boundaries/eth_admin1.shp
#
# Outputs (K:/ETH TOWERS/results/):
#   11_equity_priority_score.tif        continuous equity-adjusted score
#   11_equity_priority_classes.tif      4-class equity priority
#   11_regional_deployment_debt.csv     per-region tabular summary
#   figures/fig12_equity_priority_map.png
#   figures/fig13_regional_debt_ranking.png
#   figures/fig14_equity_vs_priority_scatter.png
#
# Runtime: ~5 minutes.
# ============================================================================

required <- c("terra", "sf", "dplyr", "tibble", "readr", "tidyr", "fs",
              "ggplot2", "patchwork", "tidyterra", "viridis", "scales",
              "stringr")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tibble)
  library(readr); library(tidyr); library(fs)
  library(ggplot2); library(patchwork); library(tidyterra)
  library(viridis); library(scales); library(stringr)
})

sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.6)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     <- "K:/ETH TOWERS"
GEE_DIR  <- file.path(ROOT, "ETH_towers")
RES_DIR  <- file.path(ROOT, "results")
FIG_DIR  <- file.path(RES_DIR, "figures")
ADM_DIR  <- file.path(ROOT, "eth_admin_boundaries")
dir_create(FIG_DIR)

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
save_fig <- function(p, file, w, h) ggsave(file, p, width = w, height = h,
                                           dpi = 300, bg = "white")

load1 <- function(p) { r <- rast(p); if (nlyr(r) > 1) r <- r[[1]]; r }
single <- function(r) if (nlyr(r) > 1) r[[1]] else r

# ---------------------------------------------------------------------------
# 1. Load all rasters and align to suitability grid
# ---------------------------------------------------------------------------
cat("=== 1. Load rasters ===\n")
suit       <- load1(file.path(RES_DIR, "07_suitability_1km.tif"))
prio_score <- load1(file.path(RES_DIR, "08_priority_score.tif"))
green      <- load1(file.path(RES_DIR, "08_priority_greenfield.tif"))
densif     <- load1(file.path(RES_DIR, "08_priority_densification.tif"))
internet   <- load1(file.path(GEE_DIR, "ETH_dhs_internet_users.tif"))
elec       <- load1(file.path(GEE_DIR, "ETH_dhs_electricity_access.tif"))
pop        <- load1(file.path(GEE_DIR, "ETH_population.tif"))

cat("Aligning all layers to 1 km grid ...\n")
internet <- resample(internet, suit, method = "bilinear", threads = TRUE)
elec     <- resample(elec,     suit, method = "bilinear", threads = TRUE)
pop      <- resample(pop,      suit, method = "bilinear", threads = TRUE)
green    <- resample(green,    suit, method = "near",     threads = TRUE)
densif   <- resample(densif,   suit, method = "near",     threads = TRUE)
prio_score <- resample(prio_score, suit, method = "bilinear", threads = TRUE)

# ---------------------------------------------------------------------------
# 2. Build digital-equity-adjusted priority surface
# ---------------------------------------------------------------------------
cat("\n=== 2. Build digital-equity-adjusted priority ===\n")

rescale01 <- function(r) {
  r <- single(r)
  v <- values(r, na.rm = TRUE)
  lo <- as.numeric(quantile(v, 0.01, na.rm = TRUE))
  hi <- as.numeric(quantile(v, 0.99, na.rm = TRUE))
  rr <- (r - lo) / (hi - lo)
  clamp(rr, 0, 1)
}

# Equity gap weights: 1 - normalised access (so high gap = high weight)
internet_gap <- 1 - rescale01(internet)
elec_gap     <- 1 - rescale01(elec)

# Combined equity weight: average of internet gap and electricity gap.
# Each gap captures a different vector of digital exclusion - mobile/internet
# uptake reflects realised connectivity, while electricity access reflects
# the structural ability to charge devices and run base-station equipment.
equity_weight <- (internet_gap + elec_gap) / 2
names(equity_weight) <- "equity_weight"

# Equity-adjusted priority = priority_score * equity_weight, then rescaled
prio_eq_raw <- single(rescale01(prio_score) * equity_weight)
prio_eq     <- rescale01(prio_eq_raw)
names(prio_eq) <- "equity_priority"

writeRaster(prio_eq, file.path(RES_DIR, "11_equity_priority_score.tif"),
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "TILED=YES", "PREDICTOR=2"))

# 4-class equity priority via quantiles
qs <- as.numeric(quantile(values(prio_eq, na.rm = TRUE),
                          probs = c(0.50, 0.75, 0.90, 1.00),
                          na.rm = TRUE))
prio_eq_cls <- classify(prio_eq,
                        rcl = matrix(c(-Inf,    qs[1], 0,
                                       qs[1],   qs[2], 1,
                                       qs[2],   qs[3], 2,
                                       qs[3],   Inf,   3),
                                     ncol = 3, byrow = TRUE),
                        include.lowest = TRUE)
prio_eq_cls <- as.factor(prio_eq_cls)
levels(prio_eq_cls) <- data.frame(value = 0:3,
                                  label = c("Below median",
                                            "Above median",
                                            "High",
                                            "Top decile"))
names(prio_eq_cls) <- "equity_class"
writeRaster(prio_eq_cls, file.path(RES_DIR, "11_equity_priority_classes.tif"),
            overwrite = TRUE, datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# Plot equity-adjusted priority
p_eq <- ggplot() +
  geom_spatraster(data = prio_eq) +
  scale_fill_viridis_c(name = "Equity-adjusted\npriority",
                       option = "magma", na.value = "transparent",
                       limits = c(0, 1)) +
  labs(title = "Digital-equity-adjusted deployment priority, Ethiopia",
       subtitle = "Modelled priority weighted by DHS-derived under-access (internet + electricity)") +
  coord_sf() + theme_pub() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 8))
save_fig(p_eq, file.path(FIG_DIR, "fig12_equity_priority_map.png"), 8, 7)

# Equity-vs-priority scatter to demonstrate orthogonality
set.seed(42)
scat <- as.data.frame(c(prio_score, equity_weight, pop), na.rm = TRUE) %>%
  as_tibble()
names(scat) <- c("priority", "equity", "pop")
scat <- scat %>% slice_sample(n = min(40000, nrow(scat)))

p_scat <- ggplot(scat, aes(priority, equity)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       labels = label_number(accuracy = 1)) +
  geom_smooth(method = "loess", se = FALSE, colour = "white",
              linewidth = 0.7) +
  labs(title = "Modelled siting priority is largely orthogonal to digital-equity gap",
       subtitle = "Hex-binned 1 km pixels (n shown by colour). High siting priority does not automatically address digital exclusion.",
       x = "Modelled priority score (suitability \u00d7 demand \u00d7 coverage gap)",
       y = "Digital-equity gap (1 = high under-access)",
       fill = "Pixels (log)") +
  theme_pub()
save_fig(p_scat, file.path(FIG_DIR, "fig14_equity_vs_priority_scatter.png"),
         8, 5.5)

# ---------------------------------------------------------------------------
# 3. Regional deployment debt
# ---------------------------------------------------------------------------
cat("\n=== 3. Regional deployment debt by ADM1 ===\n")

adm1 <- st_read(file.path(ADM_DIR, "eth_admin1.shp"), quiet = TRUE)

# HDX shapefiles use slightly different column names depending on year.
# Pick the most likely ADM1 name column.
name_candidates <- c("ADM1_EN", "ADM1_NAME", "REGIONNAME",
                     "REGION", "NAME_1", "shapeName")
adm1_namecol <- intersect(name_candidates, names(adm1))[1]
if (is.na(adm1_namecol)) {
  message("Couldn't identify ADM1 name column. Available columns:")
  print(names(adm1))
  adm1$REGION <- paste0("Region_", seq_len(nrow(adm1)))
  adm1_namecol <- "REGION"
}
adm1$REGION <- adm1[[adm1_namecol]]
cat("Using ADM1 name column: ", adm1_namecol, "\n", sep = "")

adm1_v <- vect(st_transform(adm1, crs(suit)))

# Per-region zonal statistics
extract_zonal <- function(r, fun) {
  z <- terra::extract(r, adm1_v, fun = fun, na.rm = TRUE, ID = TRUE,
                      raw = FALSE)
  z[, 2]
}

cat("Computing zonal statistics ...\n")

# Pixel area in km^2 (1km grid -> ~1 km^2 per cell after projection;
# compute exactly via cellSize)
cell_area_km2 <- cellSize(suit, unit = "km")
n_pix <- extract_zonal(!is.na(suit) * 1, "sum")          # cell count
area_km2 <- extract_zonal(cell_area_km2, "sum")
green_n  <- extract_zonal(green,  "sum")
densif_n <- extract_zonal(densif, "sum")
mean_suit <- extract_zonal(suit,  "mean")
mean_eq   <- extract_zonal(equity_weight, "mean")
mean_eq_priority <- extract_zonal(prio_eq, "mean")
total_pop <- extract_zonal(pop,   "sum")

# Population residing in greenfield-priority pixels
pop_in_green  <- extract_zonal(pop * green,  "sum")
pop_in_densif <- extract_zonal(pop * densif, "sum")

debt <- tibble(
  region              = adm1$REGION,
  area_km2            = round(area_km2, 0),
  total_pop           = round(total_pop, 0),
  greenfield_pixels   = green_n,
  densification_pixels= densif_n,
  pop_greenfield      = round(pop_in_green, 0),
  pop_densif          = round(pop_in_densif, 0),
  pct_pop_greenfield  = round(100 * pop_in_green / pmax(total_pop, 1), 2),
  pct_pop_densif      = round(100 * pop_in_densif / pmax(total_pop, 1), 2),
  greenfield_per_kkm2 = round(1000 * green_n  / pmax(area_km2, 1), 2),
  densif_per_kkm2     = round(1000 * densif_n / pmax(area_km2, 1), 2),
  mean_suitability    = round(mean_suit, 3),
  mean_equity_gap     = round(mean_eq, 3),
  mean_equity_priority= round(mean_eq_priority, 3)
) %>%
  arrange(desc(pop_greenfield))

write_csv(debt, file.path(RES_DIR, "11_regional_deployment_debt.csv"))
cat("Top 5 regions by population in greenfield-priority pixels:\n")
print(debt %>% select(region, total_pop, pop_greenfield, pct_pop_greenfield,
                      greenfield_per_kkm2, mean_equity_priority) %>%
      slice_head(n = 5))

# ---------------------------------------------------------------------------
# 4. Regional ranking figure
# ---------------------------------------------------------------------------
cat("\n=== 4. Regional ranking figure ===\n")

# Top regions ranked by greenfield population - the headline operational
# quantity (where the most underserved population could gain coverage).
debt_plot <- debt %>%
  filter(total_pop > 50000) %>%        # exclude small/empty admin units
  slice_max(pop_greenfield, n = 11) %>%
  mutate(region = factor(region, levels = rev(region)))

p_debt <- ggplot(debt_plot,
                 aes(x = pop_greenfield / 1e6, y = region,
                     fill = mean_equity_priority)) +
  geom_col(width = 0.7, alpha = 0.95) +
  geom_text(aes(label = sprintf("%.1fM (%.1f%%)",
                                pop_greenfield / 1e6, pct_pop_greenfield)),
            hjust = -0.05, size = 3.2) +
  scale_fill_viridis_c(name = "Mean equity-\nadjusted priority",
                       option = "magma", limits = c(0, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Regional deployment debt: population in greenfield-priority pixels",
       subtitle = "Bars show absolute population (millions); text shows share of regional total. Fill = mean equity-adjusted priority.",
       x = "Population in greenfield-priority pixels (millions)",
       y = NULL) +
  theme_pub()
save_fig(p_debt, file.path(FIG_DIR, "fig13_regional_debt_ranking.png"),
         10, 7)

cat("\nAll outputs in: ", RES_DIR, "\n", sep = "")
cat("DONE.\n")
