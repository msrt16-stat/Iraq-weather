# Iraq-weather

**Reproducibility repository for:** *Interpolation Gains, Extrapolation Losses: Spatial versus Global Machine Learning Models for Predicting Humidity Across Iraq's Governorates* (submitted to *Computers & Geosciences*).

This repository contains the data, R code and archived outputs needed to reproduce every model, table and figure in the paper: six models (OLS, RF, SVR and their spatial counterparts SR, SRF, SSVR), the leakage-free spatial-lag construction, random-split and leave-one-governorate-out (LOCO) validation, rank-based tests, and sensitivity analyses.

| | |
|---|---|
| Program name | Iraq-weather |
| Language | R (≥ 4.2) |
| Size | about 3,400 lines of R in 10 scripts; < 10 MB with data and outputs |
| Hardware | standard desktop/laptop; a multi-core CPU shortens scripts 2 and 8 |
| Contact | Mohsen Mohammadzadeh, Department of Statistics, Tarbiat Modares University, Tehran, Iran — mohsen_m@modares.ac.ir |
| Licence | code: MIT (see `LICENSE`); data: see [Data sources and licences](#data-sources-and-licences) |

---

## 1. Quick start

```bash
git clone https://github.com/msrt16-stat/Iraq-weather.git
cd Iraq-weather
Rscript get_boundaries.R      # downloads the GADM governorate shapefile (needed only for the maps)
Rscript run_all.R quick       # ~1-5 min: reproduces Table 2 and the variable-importance figure
Rscript run_all.R             # full pipeline: several hours
```

`run_all.R` checks that the required packages are installed, runs the scripts in the order below, writes a timing log (`outputs/run_all_log.txt`) and records the R session (`outputs/sessionInfo.txt`).

No path needs to be edited. Each script locates the project through the environment variable `IRAQ_WEATHER_ROOT`, which `run_all.R` sets to the repository root. To run a single script interactively, start R in the repository root (or call `Sys.setenv(IRAQ_WEATHER_ROOT = "/path/to/Iraq-weather")`) and `source("scripts/<name>.R")`.

### Required R packages

```r
install.packages(c("openxlsx", "readxl", "spgwr", "GWmodel", "sp", "sf", "spdep",
                   "e1071", "randomForest", "RANN", "ggplot2", "dplyr", "scatterpie",
                   "ggspatial", "ggnewscale", "gstat", "rnaturalearth",
                   "rnaturalearthdata", "blockCV"))
# for get_boundaries.R only:
install.packages(c("geodata", "terra"))
```

The package versions used for the archived outputs are listed in `outputs/sessionInfo.txt`.

---

## 2. Data

`Dalya5x.xlsx` — 2,160 rows = 18 governorates × 120 months (January 2015 – December 2024), one sheet, no missing values.

| Column | Meaning | Unit | NASA POWER parameter |
|---|---|---|---|
| `City` | governorate ID (1–18; see table below) | – | – |
| `Humidity` | specific humidity at 2 m (response) | g/kg | QV2M |
| `Temp.` | air temperature at 2 m | °C | T2M |
| `Top of Atmosphere` | top-of-atmosphere shortwave irradiance | kWh m⁻² day⁻¹ | TOA_SW_DWN |
| `Wind Speed` | wind speed | m/s | WS2M / WS10M |
| `Precipitation` | precipitation | mm/day | PRECTOTCORR |
| `Soil Moisture` | surface soil wetness | 0–1 (dimensionless) | GWETTOP |
| `Latitude`, `Longitude` | governorate coordinates (WGS84) | degrees | – |

**Row order matters.** The file has no date column. Within each `City`, rows are in chronological order (row 1 = January 2015, row 120 = December 2024). The temporal designs (rolling-origin and spatiotemporal LOCO) derive the year and month from this order. Do not sort the file.

The `City` column is stored as Excel formulas (`=A2+1`, …) with cached values; R's `openxlsx`/`readxl` read the cached values correctly.

### Governorate IDs

The ID-to-name mapping below was verified with a point-in-polygon test of each `City` coordinate against the GADM level-1 boundaries.

| City | Governorate | City | Governorate |
|---|---|---|---|
| 1 | Dihok | 10 | Wasit |
| 2 | Ninawa | 11 | Karbala |
| 3 | Arbil | 12 | Babil |
| 4 | Sulaymaniyah | 13 | Najaf |
| 5 | Kirkuk | 14 | Al-Qadisiyah |
| 6 | Salah ad-Din | 15 | Dhi Qar |
| 7 | Diyala | 16 | Maysan |
| 8 | Baghdad | 17 | Muthanna |
| 9 | Al-Anbar | 18 | Al-Basrah |

The code and CSV outputs use the raw-data column name `City`; the paper uses "governorate". They refer to the same 18 units.

---

## 3. Scripts, outputs and what they reproduce

Run in this order (`run_all.R` does this automatically). Scripts 4 and 7 read outputs of earlier scripts.

| # | Script | Main outputs (`outputs/NN_name/`) | Reproduces in the paper | Approx. runtime* |
|---|---|---|---|---|
| 1 | `weather_R_FINAL.R` | `table2_combined.csv`, `k_sensitivity_summary.csv` | Table 2; Table A1 (k-sensitivity) | minutes |
| 2 | `gwr_blockcv_scenarios.R` | `loco_spatial_by_city.csv`, `loco_spatiotemporal_by_city.csv`, `rolling_temporal_by_year.csv` | Table 3; Tables A4, A5, A7 | hours |
| 3 | `spotcheck_bisquare_12foldCV.R` | `spotcheck_spatial_loco_bisquare.csv` | Table A3 (LOCO rows) | ~1 h |
| 4 | `rerun_sr_random_split_75_25.R` | `random_split_100reps.csv`, `spotcheck_random_bisquare.csv` | random-split rows of Table 3 and Table A3 | ~1 h |
| 5 | `srf_spatial_loco.R` | `srf_spatial_loco.csv`, `srf_random_split_baseline.csv` | Table 4 (SRF); Table A6 | ~30 min |
| 6 | `ssvr_spatial_loco.R` | `ssvr_spatial_loco.csv`, `ssvr_random_split_baseline.csv` | Table 4 (SSVR); Table 5; Table A6 | ~30 min |
| 7 | `global_models_spatial_loco.R` | `rf_*`, `svr_*`, `rf_vs_srf_loco_comparison.csv`, `svr_vs_ssvr_loco_comparison.csv` | Table 4 (RF, SVR, paired tests); Table A6 | ~30 min |
| 8 | `sensitivity_analysis.R` | `sensitivity_full_summary.csv` and per-configuration files | Table A2 | hours |
| 9 | `rebuild_srf_and_importance.R` | `variable_importance.pdf` | Fig. 4 | < 1 min |
| 10 | `maps_plotting.R` | `Iraq_Climate_Interpolated_2015.png` … `2024.png` | Figs. 1–3 | minutes |

\*Single-core desktop; indicative only. Long-running scripts are checkpointed: if interrupted, re-run them and they resume from the last completed fold or repetition.

All random steps use fixed seeds (e.g. `Seed = 123`, `set.seed(20240715 + r)`), so outputs are reproducible up to platform-level floating-point differences.

### Checking a run against the archived results

The `outputs/` folder contains the results reported in the paper. After a run, compare the new CSVs with those archived in Git, for example `git diff --stat outputs/` or, in R, `all.equal(read.csv(new), read.csv(archived), tolerance = 1e-6)`.

---

## 4. Method notes relevant to the code

- **Leakage-free spatial lag.** The train/test split is made *before* any spatial-lag feature is built. Lags are computed at governorate level from training-only governorate means and a k-nearest-neighbour structure over the 18 governorate coordinates (k = 5; a governorate is never its own neighbour). In LOCO the lags are rebuilt inside every fold from the 17 training governorates only. See `weather_R_FINAL.R` (header and the function `build_city_lag()`) and Section 2.4 of the paper.
- **Statistical tests.** One-sided Mann–Whitney U (different fold sets) and one-sided Wilcoxon signed-rank (paired governorates), with rank-biserial effect sizes and 10,000-resample bootstrap percentile CIs.

---

## 5. Repository layout

```
Iraq-weather/
├── README.md
├── LICENSE
├── CITATION.cff
├── run_all.R
├── get_boundaries.R
├── Dalya5x.xlsx
├── Iraq_shp/            (created by get_boundaries.R; not version-controlled)
├── scripts/             (10 analysis scripts, see Section 3)
└── outputs/
    ├── 01_weather_R_FINAL/ … 10_maps_plotting/
    ├── run_all_log.txt  (written by run_all.R)
    └── sessionInfo.txt  (written by run_all.R)
```

---

## 6. Data sources and licences

- **Climate data:** NASA POWER Project (https://power.larc.nasa.gov/), derived from MERRA-2. NASA POWER data are freely available; please acknowledge the NASA Langley Research Center POWER Project.
- **Boundaries:** GADM v3.6 (https://gadm.org). GADM data may be used for non-commercial purposes but may not be redistributed, so they are downloaded by `get_boundaries.R` rather than stored here.
- **Code:** MIT licence.

## 7. How to cite

See `CITATION.cff`. Please cite the paper when using this code.
