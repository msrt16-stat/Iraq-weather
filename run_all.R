## =====================================================================
## run_all.R -- run the complete pipeline in the order used for the paper
##
## Usage (from the repository root):
##     Rscript run_all.R            # full pipeline (several hours)
##     Rscript run_all.R quick      # fast check: scripts 1 and 9 only
##
## Every script reads the data file and writes to outputs/<NN_name>/
## relative to the repository root, which is passed to the scripts
## through the IRAQ_WEATHER_ROOT environment variable.
## =====================================================================

root <- normalizePath(getwd())
if (!file.exists(file.path(root, "Dalya5x.xlsx"))) {
  stop("Run this script from the repository root (Dalya5x.xlsx not found in ",
       root, ").")
}
Sys.setenv(IRAQ_WEATHER_ROOT = root)

## ---- package check ------------------------------------------------------
pkgs <- c("openxlsx", "readxl", "spgwr", "GWmodel", "sp", "sf", "spdep",
          "e1071", "randomForest", "RANN", "ggplot2", "dplyr", "scatterpie",
          "ggspatial", "ggnewscale", "gstat", "rnaturalearth",
          "rnaturalearthdata", "blockCV")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       "\nInstall them with install.packages(c(\"",
       paste(missing, collapse = "\", \""), "\"))")
}

scripts <- c(
  "weather_R_FINAL.R",              # 1  Table 2, k-sensitivity
  "gwr_blockcv_scenarios.R",        # 2  SR robustness (LOCO, rolling-origin)
  "spotcheck_bisquare_12foldCV.R",  # 3  SR bi-square spot-check
  "rerun_sr_random_split_75_25.R",  # 4  SR random-split baselines
  "srf_spatial_loco.R",             # 5  SRF LOCO
  "ssvr_spatial_loco.R",            # 6  SSVR LOCO
  "global_models_spatial_loco.R",   # 7  RF/SVR LOCO + paired tests
  "sensitivity_analysis.R",         # 8  seed / grid / IDW sensitivity
  "rebuild_srf_and_importance.R",   # 9  variable-importance figure
  "maps_plotting.R"                 # 10 annual humidity maps
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] == "quick") scripts <- scripts[c(1, 9)]

log_file <- file.path(root, "outputs", "run_all_log.txt")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
cat("Run started:", format(Sys.time()), "\n", file = log_file)

for (s in scripts) {
  message("\n==== Running ", s, " ====")
  t0 <- Sys.time()
  setwd(root)                                  # scripts change the working directory
  source(file.path(root, "scripts", s), echo = FALSE, local = new.env())
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat(sprintf("%-34s finished in %6.1f min\n", s, dt), file = log_file, append = TRUE)
}
setwd(root)

## Record the computational environment used for this run
writeLines(capture.output(sessionInfo()), file.path(root, "outputs", "sessionInfo.txt"))
message("\nDone. Timing log: outputs/run_all_log.txt; environment: outputs/sessionInfo.txt")
