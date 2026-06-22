# ============================================================================
# utils.R - shared helpers used across the pipeline
# ----------------------------------------------------------------------------
# Source this near the top of any pipeline script that needs the helpers.
# ============================================================================

# ---------- Rescale a SpatRaster to [0, 1] using configurable percentiles ----------
rescale01 <- function(r, lo = 0.01, hi = 0.99) {
  v <- terra::values(r, na.rm = TRUE)
  a <- as.numeric(stats::quantile(v, lo))
  b <- as.numeric(stats::quantile(v, hi))
  rr <- (r - a) / (b - a)
  terra::clamp(rr, 0, 1)
}

# ---------- Load a single-band raster, drop extras if multi-band ----------
load_single <- function(path) {
  r <- terra::rast(path)
  if (terra::nlyr(r) > 1) r <- r[[1]]
  r
}

# ---------- Align a raster to a reference grid (CRS + extent + resolution) ----------
load_aligned <- function(path, ref) {
  r <- load_single(path)
  if (!terra::compareGeom(r, ref, stopOnError = FALSE)) {
    r <- terra::resample(r, ref, method = "bilinear", threads = TRUE)
  }
  r
}

# ---------- Find a file that may live in several known locations ----------
find_file <- function(name, candidates, optional = FALSE) {
  for (p in candidates) {
    if (file.exists(p)) {
      message(sprintf("  [found] %-25s -> %s", name, basename(p)))
      return(p)
    }
  }
  if (optional) {
    message(sprintf("  [skip ] %-25s (optional, not found)", name))
    return(NULL)
  }
  stop(sprintf("Required file '%s' not found in any candidate path.", name))
}

# ---------- Publication-quality ggplot theme ----------
theme_pub <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92",
                                               linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# ---------- Save a figure at 300 dpi with a transparent-safe background ----------
save_fig <- function(p, file, w, h, dpi = 300) {
  ggplot2::ggsave(file, p, width = w, height = h, dpi = dpi, bg = "white")
}

# ---------- Checkpoint pattern: skip if already computed ----------
checkpoint <- function(name, fn, dir = "results/checkpoints") {
  fs::dir_create(dir)
  f <- file.path(dir, paste0(name, ".rds"))
  if (file.exists(f)) {
    message(sprintf("  [resume] loading %s from checkpoint", name))
    return(readRDS(f))
  }
  result <- fn()
  saveRDS(result, f)
  gc(verbose = FALSE)
  result
}
