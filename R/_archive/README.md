# Archive: superseded script versions

These files are previous iterations of the pipeline scripts. They are kept here for full reproducibility transparency — they document the methodological decisions made during the project and the bugs that were caught in earlier runs.

**For reproducing the published results, use only the scripts in the parent `R/` directory.** The archive is read-only history.

## Pseudo-absence generation

| File | Status | Why superseded |
|---|---|---|
| `01b_regenerate_pseudo_absences.R` | Archived | Earlier rebuild after discovering ~52% of original pseudo-absences fell outside the GEE predictor extent (the bounding box covered parts of Somalia, Eritrea, Sudan, Kenya, Djibouti, and the Gulf of Aden). This version constrained sampling to a proper Ethiopian mask. |

**Current version:** `R/02_pseudo_absences.R` — adds population-weighted sampling (target-group background, Phillips et al. 2009) on top of the constrained mask, which forces the classifier to discriminate among similar populated locations rather than between populated and unpopulated cells.

## Training table assembly

| File | Status | Why superseded |
|---|---|---|
| `02_build_training_table_v2.R` | Archived | Used bulk multi-layer raster extraction which can return a different row count than the input points due to an edge case in `terra` bilinear interpolation. |

**Current version:** `R/03_build_training_table.R` (v3) — extracts predictors one layer at a time with explicit `nrow()` checks at every step, so any mismatch fails loudly rather than silently.

## ML pipeline iterations

| File | Status | Why superseded |
|---|---|---|
| `03_run_ml_pipeline_v1.R` | Archived | Earliest version, used Ookla and water_occurrence predictors that turned out to be unusable (95% NA and all-zero, respectively). |
| `03_run_ml_pipeline_v2_clean.R` | Archived | Excluded unusable predictors, but still included `ETH_dist_to_existing_tower` which caused textbook target leakage (the 1 km pseudo-absence buffer made this feature trivially separable, inflating AUC to 1.0). |
| `03_run_ml_pipeline_v3.R` | Archived | Removed `ETH_dist_to_existing_tower` from training. Stacked-ensemble metrics were still suspiciously close to 1.0 because the stack was fit on the full training set and then evaluated per fold — the meta-learner had effectively seen the test fold during its earlier training. |
| `03_run_ml_pipeline_v4.R` | Archived | Introduced honest nested out-of-fold stacking. National prediction at 100 m was slow (~7 hours), so it was sped up to 1 km in v4_final. |
| `v4_section13_14_fast.R` | Archived | Drop-in patch that sped up the national prediction step; subsequently merged directly into v4_final. |

**Current version:** `R/04_run_ml_pipeline.R` (v4_final) — all fixes baked in:
1. Honest nested OOF stacking for spatial CV metrics
2. `ETH_dist_to_existing_tower` excluded from training (target leak)
3. `ETH_ookla` and `ETH_water_occurrence` excluded (no signal / 95% NA)
4. Fast national prediction at 1 km using LightGBM (~20 min)

## Reviewer-response sensitivity suite

| File | Status | Why superseded |
|---|---|---|
| `10_reviewer_response_v1.R` | Archived | First attempt. Section 1 expected `training_table.csv` at the wrong path and Section 4 expected `08_priority_score_1km.tif` which was actually named `08_priority_score.tif`. |
| `10_patch_fix_S1_S4.R` | Archived | Hot-fix patch that added auto-detection for the two missing files. Functionality folded into v2. |
| `10b_section4_fix.R` | Archived | Standalone re-run of Section 4 with auto-detected priority filename. Folded into v2. |

**Current version:** `R/08_reviewer_response.R` (v2) — single self-contained script with all auto-detection and RAM optimizations. Addresses reviewer comments #4, #5, #6.1, #8, #19.
