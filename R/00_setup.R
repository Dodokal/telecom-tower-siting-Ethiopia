# ============================================================================
# 00_setup.R - environment bootstrap
# ----------------------------------------------------------------------------
# Run this once before anything else. It:
#   1. Checks R version
#   2. Installs all required CRAN packages (skips any already installed)
#   3. Loads them
#   4. Writes sessionInfo() to data/sessionInfo.txt for reproducibility
# ----------------------------------------------------------------------------

# ---------- 1. R version check ----------
r_ok <- getRversion() >= "4.3.0"
if (!r_ok) {
  warning("R version is below 4.3. The pipeline was developed on R 4.4.1.\n",
          "Things may still work but no guarantees.")
} else {
  message(sprintf("R version OK: %s", getRversion()))
}

# ---------- 2. Required packages ----------
required <- c(
  # Spatial core
  "terra", "sf", "raster", "stars",
  # Spatial CV and diagnostics
  "blockCV", "spdep", "gstat", "automap",
  # ML
  "tidymodels", "ranger", "xgboost", "lightgbm", "bonsai",
  "maxnet", "stacks", "glmnet",
  # SHAP / interpretation
  "fastshap", "pdp", "vip",
  # Data wrangling
  "dplyr", "tidyr", "tibble", "readr", "purrr", "fs", "yaml",
  # Visualisation
  "ggplot2", "patchwork", "viridis", "scales", "ggspatial",
  # Metrics
  "yardstick", "pROC",
  # Utilities
  "data.table"
)

to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install) > 0) {
  message("Installing missing packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# ---------- 3. Load them all ----------
suppressPackageStartupMessages({
  invisible(lapply(required, library, character.only = TRUE))
})
message("All packages loaded.")

# ---------- 4. terra and sf options ----------
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 0, memfrac = 0.5)
data.table::setDTthreads(2)

# ---------- 5. Project paths ----------
PROJECT_ROOT <- if (interactive()) getwd() else dirname(parent.frame(2)$ofile %||% ".")
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
RESULTS_DIR  <- file.path(PROJECT_ROOT, "results")
FIG_DIR      <- file.path(PROJECT_ROOT, "figures")
fs::dir_create(c(DATA_DIR, RESULTS_DIR, FIG_DIR))

# ---------- 6. Load config files ----------
HP    <- yaml::read_yaml(file.path(DATA_DIR, "hyperparameters.yaml"))
SEEDS <- yaml::read_yaml(file.path(DATA_DIR, "seeds.yaml"))
message(sprintf("Loaded %d hyperparameter blocks and %d seeds",
                length(HP), length(SEEDS)))

# ---------- 7. Write sessionInfo for reproducibility ----------
writeLines(capture.output(sessionInfo()),
           file.path(DATA_DIR, "sessionInfo.txt"))
message("sessionInfo written to data/sessionInfo.txt")

message("\n  Environment ready. Now run R/01_assemble_presence.R\n")
