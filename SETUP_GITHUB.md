# GitHub setup — step by step

Follow these in order. Open **Git Bash** on your Windows laptop and copy-paste each command. Everything between fenced code blocks is a literal command. Lines starting with `#` are comments.

---

## Step 1. Verify Git is installed

```bash
git --version
```

You should see `git version 2.x.x`. If not, install from https://git-scm.com/download/win first.

---

## Step 2. Configure Git (one-time only)

If you've never used Git on this laptop before:

```bash
git config --global user.name "Kalid Hassen Yasin"
git config --global user.email "kalid.yasin@plus.ac.at"
git config --global init.defaultBranch main
git config --global core.autocrlf true
```

Verify:

```bash
git config --global --list
```

---

## Step 3. Create the repository on GitHub.com

1. Open https://github.com/new in your browser
2. Sign in as **Dodokal**
3. Fill in:
   - **Repository name:** `telecom-tower-siting-ethiopia`
   - **Description:** `Reproduction code for: A National-Scale Ensemble Machine Learning Framework for Telecommunication Tower Site Selection in Ethiopia`
   - **Visibility:** **Private** (you can flip it to Public on acceptance)
   - **Do NOT** tick "Add a README", "Add .gitignore", or "Choose a license" — we already have those locally
4. Click **Create repository**
5. Copy the HTTPS URL it shows you. It will look like:
   ```
   https://github.com/Dodokal/telecom-tower-siting-ethiopia.git
   ```

---

## Step 4. Download the repo skeleton I built for you

I packaged the entire skeleton as a ZIP. Download `telecom-tower-siting-ethiopia.zip`, extract it to `K:/`, so you have:

```
K:/telecom-tower-siting-ethiopia/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── R/
├── data/
├── results/
├── figures/
└── docs/
```

---

## Step 5. Initialize the local repository and push

In **Git Bash**, navigate to the repo folder:

```bash
cd "/k/telecom-tower-siting-ethiopia"
```

Initialize Git locally:

```bash
git init
```

Check what Git sees:

```bash
git status
```

You should see all the files listed in red as "untracked".

Stage everything:

```bash
git add .
```

Verify (the `.gitignore` will have excluded large files automatically):

```bash
git status
```

Make the first commit:

```bash
git commit -m "Initial commit: pipeline skeleton, configs, sensitivity results"
```

Connect your local repo to the GitHub one you created in Step 3:

```bash
git remote add origin https://github.com/Dodokal/telecom-tower-siting-ethiopia.git
```

Push to GitHub:

```bash
git branch -M main
git push -u origin main
```

When it asks for credentials:
- **Username:** Dodokal
- **Password:** use a **personal access token**, NOT your GitHub password (GitHub stopped accepting passwords for git operations in 2021).

To create a token: https://github.com/settings/tokens/new
- Note: "git push from laptop"
- Expiration: 90 days
- Scopes: tick `repo`
- Click **Generate token** and **copy it now** (you only see it once)
- Paste it as the password when Git prompts

---

## Step 6. (Already done!)

Good news bro — your actual analysis scripts are already in the `R/` folder. I sorted all 20 of your uploaded versions into the latest version of each script, with older versions kept in `R/_archive/` for transparency. No need to copy anything in.

The scripts in `R/` are:

| Script | Purpose |
|---|---|
| `00_setup.R` | Environment bootstrap |
| `00a_gee_export_rasters.js` | GEE: export predictor rasters to Drive |
| `00b_gee_predictor_extraction.js` | GEE: extract values at training points |
| `01_prepare_predictors.R` | Non-GEE layers: roads, grid, DHS, flood |
| `02_pseudo_absences.R` | Population-weighted target-group background |
| `03_build_training_table.R` | Extract all predictors at presence/absence points |
| `04_run_ml_pipeline.R` | VIF → CV → RF/XGB/LGB/MaxEnt → stacked ensemble |
| `05_rf_xgb_rasters.R` | National prediction rasters for RF and XGBoost |
| `06_equity_and_regional_debt.R` | DHS-weighted priority + regional aggregation |
| `07_fix_figures.R` | Final figure cleanup |
| `08_reviewer_response.R` | All 5 sensitivity analyses for reviewer #2 |
| `utils.R` | Shared helpers |

If you make changes to any script later:

```bash
cd "/k/telecom-tower-siting-ethiopia"
git add R/
git commit -m "Update pipeline scripts"
git push
```

---

## Step 7. Verify on GitHub

Open https://github.com/Dodokal/telecom-tower-siting-ethiopia in your browser. You should see:

- The README rendered on the front page with the badges and headline-results table
- The `R/`, `data/`, `results/`, `figures/` folders
- The CSVs and PNGs visible inside `results/` and `figures/`

---

## Step 8. Future updates

Whenever you change anything:

```bash
cd "/k/telecom-tower-siting-ethiopia"
git add .
git commit -m "Describe what changed in one short sentence"
git push
```

---

## Step 9. When the paper is accepted

1. Add the final accepted PDF and supplementary to `docs/`
2. Update the citation in README and CITATION.cff with the DOI and journal name
3. Flip the repo to **Public**:
   - Settings → General → Danger Zone → Change repository visibility → Make public
4. Add the GitHub URL to the manuscript's **Data Availability Statement**:
   > Code and reproduction materials are available at https://github.com/Dodokal/telecom-tower-siting-ethiopia under CC BY 4.0.

---

## Common gotchas

**"fatal: not a git repository"** — you're not inside the project folder. Run `cd "/k/telecom-tower-siting-ethiopia"` first.

**"Authentication failed"** — you used your GitHub password instead of a personal access token. Re-do step 5 with a token.

**"large file ... exceeds GitHub's file size limit"** — a raster slipped past `.gitignore`. Run:
```bash
git rm --cached path/to/the/file.tif
git commit -m "Remove accidentally tracked large file"
```
Then add the pattern to `.gitignore` if it's not already there, commit, and push.

**Want to delete a file already pushed?**
```bash
git rm path/to/file
git commit -m "Remove unused file"
git push
```

**Want to undo your last commit but keep changes?**
```bash
git reset --soft HEAD~1
```

---

## Repository structure when done

```
telecom-tower-siting-ethiopia/
├── README.md                       ← rendered on GitHub front page
├── LICENSE                         ← CC BY 4.0
├── CITATION.cff                    ← machine-readable citation
├── .gitignore                      ← excludes large rasters
│
├── R/                              ← 10 numbered pipeline scripts + utils.R
├── data/                           ← config + metadata (no raw rasters)
├── results/                        ← the 5 sensitivity-analysis CSVs
├── figures/                        ← Figures 3, S3, S4
└── docs/                           ← manuscript PDF on acceptance
```

When this is up and live, the reviewer-comment #12 (reproducibility) is fully closed.
