# ============================================================================
# Ethiopia Telecom Tower Siting - Fix figures 13 and 14
# ----------------------------------------------------------------------------
# Two fixes:
#   * fig13: switch x-axis from absolute population (which renders poorly
#     because WorldPop sums collapse small-town pixels) to PERCENTAGE of
#     regional population in greenfield-priority pixels. Show all 11 main
#     regions instead of filtering by absolute count. Cleaner, more
#     interpretable, and emphasises the equity story.
#   * fig14: fix the broken y-axis caused by negative outliers in the
#     equity_weight rescaling. Resample the equity_weight raster from
#     scratch using a more robust 0-1 rescaling that hard-clamps to [0,1]
#     before plotting.
#
# Inputs:
#   K:/ETH TOWERS/results/11_regional_deployment_debt.csv      (already exists)
#   K:/ETH TOWERS/results/08_priority_score.tif                (already exists)
#   K:/ETH TOWERS/ETH_towers/ETH_dhs_internet_users.tif        (already exists)
#   K:/ETH TOWERS/ETH_towers/ETH_dhs_electricity_access.tif    (already exists)
#
# Outputs (overwrites previous):
#   K:/ETH TOWERS/results/figures/fig13_regional_debt_ranking.png
#   K:/ETH TOWERS/results/figures/fig14_equity_vs_priority_scatter.png
#
# Runtime: ~2 minutes.
# ============================================================================

required <- c("terra", "dplyr", "tibble", "readr", "fs",
              "ggplot2", "viridis", "scales")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(tibble); library(readr)
  library(fs); library(ggplot2); library(viridis); library(scales)
})

terra::terraOptions(progress = 0, memfrac = 0.6)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     <- "K:/ETH TOWERS"
GEE_DIR  <- file.path(ROOT, "ETH_towers")
RES_DIR  <- file.path(ROOT, "results")
FIG_DIR  <- file.path(RES_DIR, "figures")
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
# 1. FIG 13 - Regional ranking by % of regional population in greenfield
# ---------------------------------------------------------------------------
cat("=== 1. Fig 13 (regional ranking, percentage-based) ===\n")

debt <- read_csv(file.path(RES_DIR, "11_regional_deployment_debt.csv"),
                 show_col_types = FALSE)
cat("Loaded ", nrow(debt), " regions\n", sep = "")

# Drop tiny administrative oddities (Contested, Harari which is 1 city, very
# small populations) but keep the 11 main regional states + Addis + Dire Dawa.
# Filter by greenfield_pixels > 100 instead of total_pop, which is more
# reliable.
debt_plot <- debt %>%
  filter(greenfield_pixels >= 100) %>%
  arrange(desc(pct_pop_greenfield)) %>%
  mutate(region = factor(region, levels = rev(region)))

cat("Regions retained for figure: ", nrow(debt_plot), "\n", sep = "")

p_debt <- ggplot(debt_plot,
                 aes(x = pct_pop_greenfield, y = region,
                     fill = mean_equity_priority)) +
  geom_col(width = 0.7, alpha = 0.95) +
  geom_text(aes(label = sprintf("%.1f%%", pct_pop_greenfield)),
            hjust = -0.15, size = 3.4, colour = "grey20") +
  scale_fill_viridis_c(name = "Mean equity-\nadjusted priority",
                       option = "magma",
                       limits = c(0, max(debt_plot$mean_equity_priority,
                                         na.rm = TRUE)),
                       labels = label_number(accuracy = 0.01)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20)),
                     labels = label_percent(scale = 1, accuracy = 1)) +
  labs(title = "Regional deployment debt: share of regional population in greenfield-priority pixels",
       subtitle = "Bars show the percentage of each regional state's population located in pixels flagged as greenfield-priority by the LightGBM-derived deployment surface.\nFill colour shows the equity-adjusted priority intensity (priority weighted by DHS-derived under-access).",
       x = "Population in greenfield-priority pixels (% of regional total)",
       y = NULL,
       caption = "Source: own analysis. Greenfield priority = top decile of multiplicative priority score AND >5 km from existing tower.") +
  theme_pub() +
  theme(plot.caption = element_text(size = 8, colour = "grey40"),
        plot.subtitle = element_text(size = 9, colour = "grey30"))

save_fig(p_debt, file.path(FIG_DIR, "fig13_regional_debt_ranking.png"),
         11, 7)
cat("Saved -> ", file.path(FIG_DIR, "fig13_regional_debt_ranking.png"),
    "\n", sep = "")

# Print the ranked table for the manuscript
cat("\nRanked regions by % population in greenfield-priority pixels:\n")
debt_summary <- debt_plot %>%
  arrange(desc(pct_pop_greenfield)) %>%
  select(region, pct_pop_greenfield, pct_pop_densif,
         mean_suitability, mean_equity_gap, mean_equity_priority)
print(debt_summary)

# ---------------------------------------------------------------------------
# 2. FIG 14 - Equity-vs-priority scatter (with proper rescale)
# ---------------------------------------------------------------------------
cat("\n=== 2. Fig 14 (equity vs priority scatter, fixed) ===\n")

prio_score <- load1(file.path(RES_DIR, "08_priority_score.tif"))
internet   <- load1(file.path(GEE_DIR, "ETH_dhs_internet_users.tif"))
elec       <- load1(file.path(GEE_DIR, "ETH_dhs_electricity_access.tif"))

cat("Aligning DHS rasters to priority grid ...\n")
internet <- resample(internet, prio_score, method = "bilinear", threads = TRUE)
elec     <- resample(elec,     prio_score, method = "bilinear", threads = TRUE)

# Robust 0-1 rescale: clamp at fixed quantiles, hard-clip to [0,1]
rescale01_robust <- function(r, lo_q = 0.02, hi_q = 0.98) {
  r <- single(r)
  v <- values(r, na.rm = TRUE)
  lo <- as.numeric(quantile(v, lo_q, na.rm = TRUE))
  hi <- as.numeric(quantile(v, hi_q, na.rm = TRUE))
  rr <- (r - lo) / (hi - lo)
  rr <- clamp(rr, 0, 1)
  rr
}

internet_norm <- rescale01_robust(internet)
elec_norm     <- rescale01_robust(elec)
internet_gap  <- 1 - internet_norm
elec_gap      <- 1 - elec_norm
equity_weight <- (internet_gap + elec_gap) / 2
equity_weight <- clamp(equity_weight, 0, 1)
names(equity_weight) <- "equity_weight"

priority_norm <- rescale01_robust(prio_score)
names(priority_norm) <- "priority"

# Build sample for plotting
set.seed(42)
df <- as.data.frame(c(priority_norm, equity_weight), na.rm = TRUE) %>%
  as_tibble()
names(df) <- c("priority", "equity")
df <- df %>% filter(!is.na(priority), !is.na(equity),
                    priority >= 0, priority <= 1,
                    equity >= 0, equity <= 1)

cat("Pixels with both priority and equity values: ", nrow(df), "\n", sep = "")

# Compute Spearman correlation for the subtitle
cor_obj <- suppressWarnings(cor.test(df$priority, df$equity,
                                     method = "spearman", exact = FALSE))
cat("Spearman rho     : ", round(cor_obj$estimate, 4), "\n", sep = "")
cat("p-value          : ", format.pval(cor_obj$p.value, digits = 3), "\n", sep = "")

# Subsample for plotting (40k is plenty for a hex plot)
df_plot <- df %>% slice_sample(n = min(40000, nrow(df)))

# Hex bins with explicit limits
p_scat <- ggplot(df_plot, aes(priority, equity)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(option = "magma", trans = "log10",
                       labels = label_number(accuracy = 1),
                       name = "Pixels (log)") +
  geom_smooth(method = "loess", se = FALSE, colour = "white",
              linewidth = 0.7, span = 0.4) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(title = "Modelled siting priority is largely orthogonal to digital-equity gap",
       subtitle = sprintf("Spearman rho = %.3f, p = %s, n = %d 1 km pixels. High siting priority does not automatically address digital exclusion.",
                          cor_obj$estimate,
                          format.pval(cor_obj$p.value, digits = 3),
                          nrow(df)),
       x = "Modelled priority score (suitability \u00d7 demand \u00d7 coverage gap)",
       y = "Digital-equity gap (1 = high under-access)",
       caption = "Source: own analysis. Equity-gap = 1 - average of (DHS internet/mobile users, DHS electricity access), each rescaled 0-1.") +
  coord_equal() +
  theme_pub() +
  theme(plot.caption = element_text(size = 8, colour = "grey40"))

save_fig(p_scat, file.path(FIG_DIR, "fig14_equity_vs_priority_scatter.png"),
         8, 7)
cat("Saved -> ", file.path(FIG_DIR, "fig14_equity_vs_priority_scatter.png"),
    "\n", sep = "")

cat("\nDONE.\n")
