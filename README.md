# Iraq-weather
# Reproducibility Repository -- Interpolation Gains, Extrapolation Losses

This repository contains the complete R pipeline used to produce every
model, table, and robustness result reported in the paper *"Interpolation
Gains, Extrapolation Losses: Spatial versus Global Machine Learning Models
for Predicting Humidity Across Iraq's Governorates."*

## Before you run anything

- Every script sets a shared project root and writes its outputs into its
  own numbered subfolder:
  ```r
  PROJECT_ROOT <- "C:/Weather Project"
  OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "0X_scriptname")
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  setwd(OUTPUT_DIR)
  ```
  Replace `PROJECT_ROOT` with your own path before running any script.
- Place `Dalya5x.xlsx` directly in `PROJECT_ROOT`; every script reads it
  automatically (no file-selection prompt), so the pipeline can run
  non-interactively. `maps_plotting.R` also needs the Iraq governorate
  shapefile at `PROJECT_ROOT/Iraq_shp/gadm36_IRQ_1.shp`.
- A few scripts read another script's output (script 7 needs scripts 5-6;
  script 4 needs scripts 2-3), via full `PROJECT_ROOT`-based paths -- run
  the scripts in the order below at least once so these files exist.
- Required R packages (install once): `openxlsx`, `readxl`, `spgwr`,
  `GWmodel`, `sp`, `sf`, `spdep`, `e1071`, `randomForest`, `RANN`,
  `ggplot2`, `dplyr`, `scatterpie`, `ggspatial`, `ggnewscale`, `gstat`,
  `rnaturalearth`, `rnaturalearthdata`, `blockCV`, `parallel` (base R).
- Several scripts are checkpointed: they write partial results as they go
  and skip completed folds/reps on a re-run. If interrupted, just re-run
  the script -- it resumes rather than restarts.
- Runtime varies enormously, from under a minute
  (`rebuild_srf_and_importance.R`) to several hours
  (`gwr_blockcv_scenarios.R`, `sensitivity_analysis.R`). See the runtime
  warning inside each script before running.

## Recommended execution order

1. `weather_R_FINAL.R`
2. `gwr_blockcv_scenarios.R`
3. `spotcheck_bisquare_12foldCV.R`
4. `rerun_sr_random_split_75_25.R`
5. `srf_spatial_loco.R`
6. `ssvr_spatial_loco.R`
7. `global_models_spatial_loco.R`
8. `sensitivity_analysis.R`
9. `rebuild_srf_and_importance.R`
10. `maps_plotting.R`

## Scripts, their outputs, and what each supports in the paper

| # | Script | Output files (in `outputs/0#_scriptname/`) | Supports |
|---|---|---|---|
| 1 | `weather_R_FINAL.R` | `k_sensitivity_summary.csv`, `train_metrics_for_table2.csv`, `test_metrics_for_table2.csv`, `table2_combined.csv` | **Table 2**; k-sensitivity table |
| 2 | `gwr_blockcv_scenarios.R` | `loco_spatial_by_city.csv`, `loco_spatiotemporal_by_city.csv`, `rolling_temporal_by_year.csv` (+ optional `loco_summary_stats.csv`, `blockcv_comparison_results_single_fold.csv`); also `random_split_100reps.csv` (superseded, see #4) | SR robustness: Rolling-origin / Spatial / Spatiotemporal LOCO rows; per-governorate LOCO appendix |
| 3 | `spotcheck_bisquare_12foldCV.R` | `spotcheck_spatial_loco_bisquare.csv`; also an internal random-split baseline (superseded, see #4) | Bi-square/12-fold-CV spot-check table |
| 4 | `rerun_sr_random_split_75_25.R` | `random_split_100reps.csv`, `spotcheck_random_bisquare.csv` (both 75/25, **authoritative**) | Random-split baseline rows of both CV tables above |
| 5 | `srf_spatial_loco.R` | `srf_spatial_loco.csv`, `srf_random_split_baseline.csv` | SRF's random-split-vs-LOCO comparison |
| 6 | `ssvr_spatial_loco.R` | `ssvr_spatial_loco.csv`, `ssvr_random_split_baseline.csv` | SSVR's random-split-vs-LOCO comparison; Maysan negative-$R^2$ finding and neighbor table |
| 7 | `global_models_spatial_loco.R` | `rf_spatial_loco.csv`, `rf_random_split_baseline.csv`, `svr_spatial_loco.csv`, `svr_random_split_baseline.csv` (required); + `rf_vs_srf_loco_comparison.csv`, `svr_vs_ssvr_loco_comparison.csv` (derivable, optional) | Paired RF-vs-SRF / SVR-vs-SSVR tests -- the paper's central finding |
| 8 | `sensitivity_analysis.R` | `sensitivity_baseline_seed123.csv`, `sensitivity_seed456.csv`, `sensitivity_seed789.csv`, `sensitivity_wide_grid.csv`, `sensitivity_idw_weights.csv`, `sensitivity_full_summary.csv` | Robustness to seed / grid / spatial-weighting (Limitations) |
| 9 | `rebuild_srf_and_importance.R` | `variable_importance.pdf` | Variable-importance figure (Results) |
| 10 | `maps_plotting.R` | 10 PNG maps (`Iraq_Climate_Interpolated_2015.png` ... `2024.png`) | Ten-year spatio-temporal humidity figures |

Output files from previous runs can be archived directly in the matching
folder above, without re-running anything. Git does not track empty
folders, so there's no need for placeholder files in any folder you leave
empty.

## Repository layout

```
C:/Weather Project/
├── Dalya5x.xlsx
├── Iraq_shp/
│   └── gadm36_IRQ_1.shp   (plus its .dbf/.shx/.prj companions)
├── scripts/
│   ├── weather_R_FINAL.R
│   ├── gwr_blockcv_scenarios.R
│   ├── spotcheck_bisquare_12foldCV.R
│   ├── rerun_sr_random_split_75_25.R
│   ├── srf_spatial_loco.R
│   ├── ssvr_spatial_loco.R
│   ├── global_models_spatial_loco.R
│   ├── sensitivity_analysis.R
│   ├── rebuild_srf_and_importance.R
│   └── maps_plotting.R
├── outputs/
│   ├── 01_weather_R_FINAL/
│   ├── 02_gwr_blockcv_scenarios/
│   ├── 03_spotcheck_bisquare_12foldCV/
│   ├── 04_rerun_sr_random_split_75_25/
│   ├── 05_srf_spatial_loco/
│   ├── 06_ssvr_spatial_loco/
│   ├── 07_global_models_spatial_loco/
│   ├── 08_sensitivity_analysis/
│   ├── 09_rebuild_srf_and_importance/
│   └── 10_maps_plotting/
└── README.md
```

This is both the local working layout and the repository layout -- running
the scripts produces exactly this structure, so nothing needs to be
manually reorganized afterward.

## A note on governorate-name correction in `maps_plotting.R`

The Iraq governorates shapefile (`gadm36_IRQ_1.shp`) spells two
governorates as `Sala ad-Din` and `Al-Muthannia`. Since this shapefile is
a separate GIS asset that has not been edited, `maps_plotting.R` corrects
these two names to `Salah al-Din` and `Muthanna` immediately after reading
the shapefile, so every downstream use -- the City-ID lookup table and the
map's text labels -- automatically shows the corrected spelling. If you
ever replace `gadm36_IRQ_1.shp` with a shapefile that already spells these
correctly, this correction step becomes a harmless no-op.

## "City" (code/CSVs) vs. "Governorate" (paper)

Every script and output CSV identifies each of Iraq's 18 administrative
units using a numeric `City` column (1-18), inherited from the source
Excel file's own column name. The paper uses the more precise term
"governorate" instead. These refer to the exact same 18 units -- `City`
was kept as-is in the code, since it's a raw-data column name depended
upon across many interdependent scripts and CSV files, and renaming it
everywhere would add risk without changing what the data represents. Use
this table to map any `City` value to its governorate name:

| City | Governorate | City | Governorate |
|---|---|---|---|
| 1 | Dihok | 10 | Wasit |
| 2 | Ninawa | 11 | Karbala |
| 3 | Arbil | 12 | Babil |
| 4 | Sulaymaniyah | 13 | Najaf |
| 5 | Kirkuk | 14 | Maysan |
| 6 | Salah al-Din | 15 | Qadisiyah |
| 7 | Diyala | 16 | Dhi Qar |
| 8 | Baghdad | 17 | Muthanna |
| 9 | Al-Anbar | 18 | Al-Basrah |
