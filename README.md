# Machine Learning and Geospatial Modeling Reveal Telecommunication Tower Suitability, Deployment Priorities, and Digital Connectivity Gaps Across Ethiopia

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![R version](https://img.shields.io/badge/R-%E2%89%A5%204.3-blue.svg)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/status-under%20review-yellow.svg)]()

This repository contains the code, predictor metadata, and reproduction instructions for the manuscript:

> **Reja, A. A. & Yasin, K. H. (2026).** *Machine Learning and Geospatial Modeling Reveal Telecommunication Tower Suitability, Deployment Priorities, and Digital Connectivity Gaps Across Ethiopia.* Manuscript under review.

The framework integrates open Earth observation, Demographic and Health Survey (DHS), and crowdsourced infrastructure data through a stacked ensemble of LightGBM, XGBoost, Random Forest, and MaxEnt classifiers under nested spatial block cross-validation. It produces a national 1 km suitability surface, an operational priority surface, and a digital-equity-adjusted priority surface that systematically reveals geographies of digital exclusion missed by conventional siting maps.

---

## Quick links

- 📄 [Manuscript and supplementary materials](./docs/)
- 🔬 [Reproduction guide](#reproduction)
- 🛠️ [Pipeline scripts](./R/)
- 📊 [Headline results](#headline-results)
- 📚 [Data sources](#data-sources)
- ✉️ [Contact](#contact)

---

## Headline results

| Metric | Value |
|---|---|
| Tree-based AUC range (10-fold spatial CV) | **0.963–0.970** |
| LightGBM Brier score | **0.001529** |
| Greenfield-priority area | **109,630 km² (9.87% of national)** |
| Densification-priority area | **6,054 km² (0.54%)** |
| Spearman ρ (priority vs equity gap) | **−0.037** (essentially orthogonal) |
| Regional rank stability across 4 classifiers | **ρ = 0.91** |
| Variogram range of residuals | **7.45 km** (≪ 50 km block size) |

See `results/` for the full set of cross-validated metrics, Moran I diagnostics, calibration scores, and threshold-sensitivity tables.

---

## Repository structure

```
telecom-tower-siting-ethiopia/
├── README.md                       ← you are here
├── LICENSE                         ← CC BY 4.0
├── CITATION.cff                    ← machine-readable citation metadata
├── .gitignore                      ← excludes large rasters, secrets
│
├── R/                              ← analysis pipeline (numbered, sequential)
│   ├── 00_setup.R                  ← install + load packages, set paths
│   ├── 00a_gee_export_rasters.js   ← Google Earth Engine: export predictor rasters
│   ├── 00b_gee_predictor_extraction.js ← GEE: extract values at training points
│   ├── 01_prepare_predictors.R     ← non-GEE layers: roads, grid, DHS, flood
│   ├── 02_pseudo_absences.R        ← population-weighted target-group background
│   ├── 03_build_training_table.R   ← extract all predictors at presence/absence pts
│   ├── 04_run_ml_pipeline.R        ← VIF → spatial CV → RF/XGB/LGB/MaxEnt → stack
│   ├── 05_rf_xgb_rasters.R         ← national prediction surfaces for RF and XGB
│   ├── 06_equity_and_regional_debt.R ← DHS-weighted priority + regional aggregation
│   ├── 07_fix_figures.R            ← regenerate figures 13–14 with corrected axes
│   ├── 08_reviewer_response.R      ← all 5 sensitivity analyses for reviewer #2
│   ├── utils.R                     ← shared helpers
│   └── _archive/                   ← earlier script versions (transparency)
│
├── data/
│   ├── README.md                   ← data sources, licences, how to obtain
│   ├── predictor_metadata.csv      ← machine-readable predictor list
│   ├── hyperparameters.yaml        ← tuned hyperparameter values for each learner
│   └── seeds.yaml                  ← random seeds for reproducibility
│
├── results/
│   ├── README.md                   ← description of each output file
│   ├── r1_moran_results.csv        ← spatial autocorrelation diagnostics
│   ├── r1_variogram_fit.csv        ← variogram range = 7.45 km
│   ├── r2_brier_scores.csv         ← calibration metric
│   ├── r4_equity_weight_sensitivity.csv
│   └── r5_threshold_sensitivity.csv
│
├── figures/                        ← published figures (PNG, 300 dpi)
│   ├── README.md
│   ├── fig_var_importance.png      ← Figure 3, four-panel SHAP / permutation
│   ├── r1_moran_variogram.png      ← Figure S3
│   └── r2_reliability_diagram.png  ← Figure S4
│
└── docs/
    ├── manuscript.pdf              ← latest accepted version (added on acceptance)
    └── supplementary.pdf
```

---

## Reproduction

### Software requirements

- **R ≥ 4.3** (tested on R 4.4.1 on Windows 11)
- ~8 GB RAM minimum (the pipeline was developed on an 8 GB laptop; see notes below)
- ~5 GB free disk for intermediate rasters

### Setup (one-time)

Clone the repository and install dependencies:

```bash
git clone https://github.com/Dodokal/telecom-tower-siting-ethiopia.git
cd telecom-tower-siting-ethiopia
```

In R or RStudio:

```r
source("R/00_setup.R")
```

This installs the required CRAN packages and verifies your environment. The full list and its pinned versions are in `data/sessionInfo.txt`.

### Run the full pipeline

The pipeline runs in two stages. First, the predictor rasters are exported from Google Earth Engine (one-time setup). Then the R pipeline runs sequentially:

**Stage 1 — Google Earth Engine (one-time, ~2 hours)**

Open the GEE Code Editor at https://code.earthengine.google.com and paste in:

1. `R/00a_gee_export_rasters.js` — exports all GEE predictor rasters to Google Drive (`ETH_towers/` folder). Click RUN on each task.
2. `R/00b_gee_predictor_extraction.js` — extracts values at training points (optional, for cross-checking).

Download the exported GeoTIFFs from Google Drive into your local `ETH_towers/` folder.

**Stage 2 — R pipeline (sequential)**

```r
source("R/00_setup.R")                       # bootstrap environment

source("R/01_prepare_predictors.R")          # roads, grid, DHS, flood (non-GEE)
source("R/02_pseudo_absences.R")             # target-group background sampling
source("R/03_build_training_table.R")        # extract all predictors at points
source("R/04_run_ml_pipeline.R")             # VIF → CV → RF/XGB/LGB/MaxEnt → stack
source("R/05_rf_xgb_rasters.R")              # national rasters for RF and XGB
source("R/06_equity_and_regional_debt.R")    # equity-weighted priority + ADM1 stats
source("R/07_fix_figures.R")                 # final figure cleanup
source("R/08_reviewer_response.R")           # 5 sensitivity analyses (reviewer #2)
```

End-to-end runtime is roughly 6–10 hours on an 8 GB laptop (most of it in `04_run_ml_pipeline.R` due to spatial-CV hyperparameter tuning).

### Tight-RAM mode (8 GB or less)

The default settings are already RAM-conservative. If you still hit OOM errors:

1. Reduce `n_threads` in `data/hyperparameters.yaml` from 2 to 1
2. Set `terra::terraOptions(memfrac = 0.3)` at the top of any heavy script
3. Skip the four-classifier stacked ensemble and use LightGBM only (the paper shows they are statistically indistinguishable; paired Wilcoxon p = 0.286). In `R/04_run_ml_pipeline.R`, set `RUN_STACK <- FALSE` near the top.

---

## Data sources

The 23 predictor variables are listed in `data/predictor_metadata.csv` with their licences and download URLs. All sources are open and reproducible:

| Group | Layer | Source |
|---|---|---|
| Demand | Population density | WorldPop 2020 constrained (CC BY 4.0) |
| Demand | Built-up surface | GHSL GHS-BUILT-S R2023A (CC BY 4.0) |
| Demand | Night-time lights | VIIRS DNB Annual V22 |
| Demand | DHS internet & electricity | DHS Program 2016 Ethiopia (registered access) |
| Terrain | Elevation, slope, aspect, TRI, TPI | SRTM v3 (public domain) |
| Access | Roads, airports, urban centres | OpenStreetMap via Geofabrik (ODbL) |
| Access | Power grid | gridfinder predictive MV grid (CC BY 4.0) |
| Environment | Land cover | ESA WorldCover 2021 v200 (CC BY 4.0) |
| Environment | Precipitation | CHIRPS v2.0 |
| Constraint | Protected areas | WDPA May 2024 |
| Constraint | Flood inundation | Global Flood Model intercomparison |

**Tower presences** were derived from OpenStreetMap features tagged `man_made = tower` with `tower:type = communication` plus `man_made = mast`. The extracted points are not redistributed here as they originate from OSM (ODbL); regenerate them with `R/01_assemble_presence.R` using the Geofabrik Ethiopia extract.

**DHS data** require free registration at https://dhsprogram.com. Place the 2016 Ethiopia cluster shapefile in `data/dhs_2016/`.

---

## Citation

If you use this code or the framework, please cite:

> Reja, A. A. & Yasin, K. H. (2026). Machine Learning and Geospatial Modeling Reveal Telecommunication Tower Suitability, Deployment Priorities, and Digital Connectivity Gaps Across Ethiopia. *[Journal]* (under review).

A machine-readable BibTeX entry is in `CITATION.cff`.

---

## Reviewer reproducibility statement

The reviewer of an earlier version of this manuscript requested specific reproducibility materials. They are all present in this repository:

| Reviewer request | Location |
|---|---|
| Code repository | This repo |
| Parameter files | `data/hyperparameters.yaml` |
| Random seeds | `data/seeds.yaml` |
| Trained models | Generated by `R/04_run_ml_pipeline.R`; saved to `results/models/` |
| Workflow diagram | Figure 2 in the main manuscript |
| Hyperparameter ranges | `data/hyperparameters.yaml` (search grids and final values) |
| Predictor preprocessing scripts | `R/01_prepare_predictors.R`, `R/03_build_training_table.R` |
| Sensitivity analyses | `R/08_reviewer_response.R` |

---

## Licence

Code is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Predictor rasters retain their original licences (see table above). Tower presence locations from OSM are redistributable under the [Open Database Licence (ODbL)](https://opendatacommons.org/licenses/odbl/).

---

## Contact

- **Kalid Hassen Yasin** (corresponding)
  Department of Geoinformatics — Z_GIS, University of Salzburg, Austria
  Department of Geography and Environmental Studies, Haramaya University, Ethiopia

- **Amira Ahmed Reja** — Department of Electrical and Computer Engineering, Adigrat University, Ethiopia

Issues and pull requests are welcome.
