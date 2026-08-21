## ============================================================
## *** ARCHIVAL NOTE (added when consolidating scripts for the GitHub
## data-availability repository) ***
##
## 1. RANDOM SPLIT RATIO: Scenario 1 below uses a 70/30 split
##    (floor(0.7 * n_total)) and is SUPERSEDED for the numbers reported
##    in the final paper. The Random-split baseline actually reported
##    (RMSE = 0.7235 +/- 0.0216) was produced by
##    rerun_sr_random_split_75_25.R using a 75/25 split instead; use
##    that script's output for the authoritative Random-split number.
##    The Rolling-origin, Spatial LOCO, and Spatiotemporal LOCO
##    scenarios below do NOT depend on this ratio and remain fully
##    valid and unchanged -- their outputs (rolling_temporal_by_year.csv,
##    loco_spatial_by_city.csv, loco_spatiotemporal_by_city.csv) are the
##    authoritative sources used throughout the paper.
##
## 2. EXTRA SCENARIOS NOT IN THE FINAL PAPER: this script also runs
##    three additional exploratory scenarios -- "Temporal block"
##    (single fixed 2015-21/2022-24 split), "Spatiotemporal block", and
##    "Spatial block only" (both via a single k-means-clustered
##    governorate hold-out) -- and saves
##    blockcv_comparison_results_single_fold.csv and
##    loco_summary_stats.csv. These reflect an earlier exploratory stage
##    of the analysis and are NOT referenced anywhere in the final
##    manuscript, which reports only the four scenarios described above
##    (per-governorate LOCO superseded the k-means cluster-block
##    approach). They are left in this archived script for transparency
##    but should not be cited as paper results.
##
## 3. PARALLEL EXECUTION: this script uses a PSOCK parallel cluster
##    (parLapply). If this causes worker crashes or connection errors on
##    a given system, gwr_blockcv_scenarios_reconstructed.R provides a
##    fully sequential, checkpointed equivalent for the same four core
##    scenarios.
## ============================================================

## ============================================================
## Block-CV sensitivity analysis for the GWR humidity model
## Scenarios: (1) Random split  (2) Temporal block  (3) Spatiotemporal block
##
## ASSUMPTION TO VERIFY BEFORE TRUSTING ANY NUMBER BELOW:
##   The file has NO Year/Month column. There are 18 cities x 120 monthly
##   rows each, and the within-city sequence shows a clear 12-month
##   seasonal cycle (Temp rises/falls repeatedly). I am assuming
##   row 1 of each city = January 2015, row 120 = December 2024
##   (10 full years). If this is wrong, every train/test boundary below
##   is wrong. Confirm against your original data source / collection log.
## ============================================================

set.seed(20240715)  # fixed seed for reproducibility, cite in the paper

## ---- 0. Packages ----
required_pkgs <- c("readxl", "dplyr", "GWmodel", "sp", "blockCV", "sf", "parallel")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

PROJECT_ROOT <- "C:/Weather Project"
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "02_gwr_blockcv_scenarios")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)

library(readxl)
library(dplyr)
library(GWmodel)
library(sp)
library(blockCV)
library(sf)
library(parallel)

rm(list=ls())

set.seed(20240715)

## ---- 1. Load & clean ----
path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
raw <- read_excel(path)  ## read the data ("Dalya5x.xlsx")

df <- raw %>%
  rename(
    City      = City,
    Humidity  = Humidity,
    Temp      = `Temp.`,
    TOA       = `Top of Atmosphere`,   # source file's header was corrected from the earlier "Top Of Amosphere" typo
    Wind      = `Wind Speed`,
    Precip    = Precipitation,
    Soil      = `Soil Moisture`,
    Lat       = Latitude,
    Lon       = Longitude
  ) %>%
  arrange(City) %>%
  group_by(City) %>%
  mutate(
    month_idx = row_number(),                 # 1..120 within each city, in file order
    Year      = 2015 + (month_idx - 1) %/% 12,   # <-- VERIFY this assumption
    Month     = ((month_idx - 1) %% 12) + 1
  ) %>%
  ungroup()

stopifnot(all(df %>% count(City) %>% pull(n) == 120))
cat("Cities:", n_distinct(df$City), " | Rows:", nrow(df),
    " | Year range assumed:", range(df$Year), "\n")

model_formula <- Humidity ~ Temp + TOA + Wind + Precip + Soil

## ---- Helper: fit GWR on train, predict on test, return metrics ----
fit_predict_gwr <- function(train_df, test_df, label) {

  train_sp <- train_df
  coordinates(train_sp) <- ~ Lon + Lat
  test_sp <- test_df
  coordinates(test_sp) <- ~ Lon + Lat

  ## Bandwidth selected ONLY on train (this is where leakage most often
  ## sneaks back in if you forget and tune on the full data)
  bw_train <- bw.gwr(
    model_formula, data = train_sp,
    approach = "AICc", kernel = "gaussian", adaptive = TRUE
  )

  pred <- gwr.predict(
    model_formula, data = train_sp,
    predictdata = test_sp,
    bw = bw_train, kernel = "gaussian", adaptive = TRUE
  )

  yhat <- pred$SDF$prediction
  yobs <- test_df$Humidity

  resid <- yobs - yhat
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  mae   <- mean(abs(resid), na.rm = TRUE)

  ## Lag-1 temporal autocorrelation of residuals, computed PER CITY then
  ## averaged -- a diagnostic for whether temporal leakage is still present
  resid_df <- test_df %>% mutate(resid = resid) %>% arrange(City, month_idx)
  lag1_by_city <- resid_df %>%
    group_by(City) %>%
    summarise(
      lag1 = if (n() > 3) cor(resid[-1], resid[-n()], use = "complete.obs") else NA_real_,
      .groups = "drop"
    )
  mean_lag1 <- mean(lag1_by_city$lag1, na.rm = TRUE)

  cat(sprintf("\n[%s] bw = %.2f | n_train = %d | n_test = %d\n",
              label, bw_train, nrow(train_df), nrow(test_df)))
  cat(sprintf("  RMSE = %.4f | MAE = %.4f | mean residual lag-1 autocorr = %.3f\n",
              rmse, mae, mean_lag1))

  list(label = label, rmse = rmse, mae = mae, lag1 = mean_lag1,
       n_train = nrow(train_df), n_test = nrow(test_df))
}

results <- list()

## ============================================================
## Parallel cluster setup (used by Scenario 1 and both LOCO loops below)
## Cross-platform: PSOCK works on Windows, macOS, and Linux alike.
## (FORK would be slightly faster on Linux/macOS but silently fails on
##  Windows, so PSOCK is the safer default for a shared script.)
## ============================================================
N_CORES <- max(1, parallel::detectCores() - 1)  # leave one core free
cat("\nUsing", N_CORES, "cores for parallel execution.\n")

cl <- makeCluster(N_CORES, type = "PSOCK", outfile = "")

## outfile = "" is important for debugging: by default, PSOCK workers send
## their stdout/stderr to /dev/null (or nowhere useful on Windows), so if a
## worker crashes while loading a package (e.g. GWmodel/sp dependency issues)
## you get an opaque "error reading/writing connection" with no real cause.
## With outfile = "", worker messages and errors print to the main R console.

## Every worker process starts empty -- it needs the packages, the data,
## the model formula, and the fit_predict_gwr() function itself.
clusterEvalQ(cl, {
  library(dplyr)
  library(GWmodel)
  library(sp)
  NULL
})
clusterExport(cl, varlist = c("df", "model_formula", "fit_predict_gwr"))

## Always stop the cluster on exit, even if something errors out later in
## the script -- otherwise orphaned R processes are left running.
on.exit(if (exists("cl")) try(stopCluster(cl), silent = TRUE), add = TRUE)

## ============================================================
## Scenario 1: RANDOM split (the ORIGINAL, flawed approach)
## Repeated over N_REPS different seeds so we get a mean +/- sd,
## comparable on equal footing to the LOCO mean +/- sd below.
## A single random split gives one number with no notion of its own
## sampling variability; that is not a fair comparison to 18-fold LOCO.
## Parallelized across repetitions with parLapply: each repetition is
## fully independent (own seed, own split, own GWR fit), so this is
## an embarrassingly-parallel loop with no shared state between workers.
## ============================================================
N_REPS <- 100
n_total <- nrow(df)  # named n_total, NOT n -- dplyr::n() would mask a bare `n` on workers
clusterExport(cl, varlist = c("n_total"))

random_rep_worker <- function(r) {
  set.seed(20240715 + r)  # distinct, reproducible seed per repetition
  rand_idx   <- sample.int(n_total, size = floor(0.7 * n_total))
  train_rand <- df[rand_idx, ]
  test_rand  <- df[-rand_idx, ]

  res <- tryCatch(
    fit_predict_gwr(train_rand, test_rand, sprintf("Random split rep=%d", r)),
    error = function(e) {
      list(label = sprintf("Random split rep=%d", r), rmse = NA_real_,
           mae = NA_real_, lag1 = NA_real_,
           n_train = nrow(train_rand), n_test = nrow(test_rand),
           error_msg = conditionMessage(e))
    }
  )
  res$rep <- r
  res
}

## RUNTIME NOTE: with N_CORES workers, 100 reps take roughly
## (100 / N_CORES) x (time per single GWR fit) instead of 100x --
## check N_CORES printed above to gauge expected speedup on your machine.
random_fold_results <- parLapply(cl, seq_len(N_REPS), random_rep_worker)

## Surface any failed reps here (parLapply doesn't message() back to the
## console live the way a sequential loop would)
failed_reps <- Filter(function(x) is.na(x$rmse), random_fold_results)
if (length(failed_reps) > 0) {
  for (fr in failed_reps) {
    message(sprintf("  [RANDOM REP FAILED] rep=%d : %s", fr$rep,
                     if (!is.null(fr$error_msg)) fr$error_msg else "unknown error"))
  }
}

random_reps_df <- bind_rows(lapply(random_fold_results, function(x) {
  data.frame(rep = x$rep, n_train = x$n_train, n_test = x$n_test,
             RMSE = x$rmse, MAE = x$mae, Lag1 = x$lag1)
}))

n_failed_random <- sum(is.na(random_reps_df$RMSE))
if (n_failed_random > 0) {
  warning(sprintf("%d of %d random-split repetitions failed to fit -- excluded from summary stats.",
                   n_failed_random, N_REPS))
}

write.csv(random_reps_df, "random_split_100reps.csv", row.names = FALSE)
cat("\nSaved: random_split_100reps.csv (per-repetition RMSE/MAE/Lag1 across", N_REPS, "seeds)\n")

valid_random <- random_reps_df %>% filter(!is.na(RMSE))

results$random <- list(
  label   = sprintf("Random split (70/30, mean of %d reps)", nrow(valid_random)),
  rmse    = mean(valid_random$RMSE),
  mae     = mean(valid_random$MAE),
  lag1    = mean(valid_random$Lag1, na.rm = TRUE),
  n_train = round(mean(valid_random$n_train)),
  n_test  = round(mean(valid_random$n_test)),
  rmse_sd = sd(valid_random$RMSE),
  mae_sd  = sd(valid_random$MAE)
)


cat(sprintf("\n[Random split, %d reps] RMSE = %.4f +/- %.4f | MAE = %.4f +/- %.4f\n",
            nrow(valid_random), results$random$rmse, results$random$rmse_sd,
            results$random$mae, results$random$mae_sd))

## ============================================================
## Scenario 2: TEMPORAL block CV
## Train: 2015-2021 (all 18 cities) | Test: 2022-2024 (all 18 cities)
## Fixes temporal leakage. Does NOT fix spatial confounding
## (same 18 locations appear in both train and test).
## Kept as-is for continuity with earlier reported numbers, but see
## Scenario 2b below for the version that should actually go in the paper.
## ============================================================
train_temporal <- df %>% filter(Year <= 2021)
test_temporal  <- df %>% filter(Year >= 2022)

results$temporal <- fit_predict_gwr(train_temporal, test_temporal, "Temporal block (2015-21 / 2022-24)")

## ============================================================
## Scenario 2b: ROLLING-ORIGIN (walk-forward) temporal CV
##
## The single fixed split above has no natural repetition axis: there is
## only one way to cut at the 2021/2022 boundary, so it cannot honestly be
## reported as mean +/- sd. The standard fix in time-series forecast
## evaluation is an expanding-window walk-forward scheme: grow the
## training window one year at a time and always test on the single next
## year. This gives multiple genuinely independent folds -- comparable in
## spirit to the 18 LOCO folds -- and a defensible mean +/- sd.
##
## Fold 1: train 2015-2019 (5y) -> test 2020
## Fold 2: train 2015-2020 (6y) -> test 2021
## Fold 3: train 2015-2021 (7y) -> test 2022   <- this is the Scenario 2 split
## Fold 4: train 2015-2022 (8y) -> test 2023
## Fold 5: train 2015-2023 (9y) -> test 2024
## ============================================================
first_test_year <- 2020   # first year used as a test fold (needs >=5 years of prior training)
last_test_year  <- 2024
test_years <- first_test_year:last_test_year

rolling_worker <- function(test_year) {
  train_fold <- df %>% filter(Year < test_year)
  test_fold  <- df %>% filter(Year == test_year)
  lbl <- paste0("Rolling-temporal test_year=", test_year)

  res <- tryCatch(
    fit_predict_gwr(train_fold, test_fold, lbl),
    error = function(e) {
      list(label = lbl, rmse = NA_real_, mae = NA_real_, lag1 = NA_real_,
           n_train = nrow(train_fold), n_test = nrow(test_fold),
           error_msg = conditionMessage(e))
    }
  )
  res$test_year <- test_year
  res
}

## Sequential version (use this if you did NOT keep the parallel cluster
## from the earlier version of the script -- this matches what you say
## ran fine in reasonable time without parallelization):
rolling_fold_results <- lapply(test_years, rolling_worker)

## If you DO have the `cl` cluster from the parallel version still set up,
## you can swap the line above for:
##   clusterExport(cl, varlist = c("df", "fit_predict_gwr"))
##   rolling_fold_results <- parLapply(cl, test_years, rolling_worker)
## but with only 5 folds the parallel overhead usually isn't worth it.

failed_rolling <- Filter(function(x) is.na(x$rmse), rolling_fold_results)
if (length(failed_rolling) > 0) {
  for (fr in failed_rolling) {
    message(sprintf("  [ROLLING FOLD FAILED] test_year=%d : %s", fr$test_year,
                     if (!is.null(fr$error_msg)) fr$error_msg else "unknown error"))
  }
}

rolling_temporal_df <- bind_rows(lapply(rolling_fold_results, function(x) {
  data.frame(test_year = x$test_year, n_train = x$n_train, n_test = x$n_test,
             RMSE = x$rmse, MAE = x$mae, Lag1 = x$lag1)
}))

cat("\n\n==================== ROLLING-ORIGIN TEMPORAL CV (5 folds) ====================\n")
print(rolling_temporal_df, row.names = FALSE)

write.csv(rolling_temporal_df, "rolling_temporal_by_year.csv", row.names = FALSE)
cat("\nSaved: rolling_temporal_by_year.csv\n")

valid_rolling <- rolling_temporal_df %>% filter(!is.na(RMSE))

results$temporal_rolling <- list(
  label   = sprintf("Temporal rolling-origin (%d folds, mean +/- sd)", nrow(valid_rolling)),
  rmse    = mean(valid_rolling$RMSE),
  mae     = mean(valid_rolling$MAE),
  lag1    = mean(valid_rolling$Lag1, na.rm = TRUE),
  n_train = round(mean(valid_rolling$n_train)),
  n_test  = round(mean(valid_rolling$n_test)),
  rmse_sd = sd(valid_rolling$RMSE),
  mae_sd  = sd(valid_rolling$MAE)
)

cat(sprintf("\n[Rolling-origin temporal, %d folds] RMSE = %.4f +/- %.4f | MAE = %.4f +/- %.4f\n",
            nrow(valid_rolling), results$temporal_rolling$rmse, results$temporal_rolling$rmse_sd,
            results$temporal_rolling$mae, results$temporal_rolling$mae_sd))

## NOTE ON SAMPLE SIZE: 5 folds is a small sample for a standard deviation --
## the sd above should be interpreted as indicative, not a tight estimate.
## Do not over-claim precision from it in the paper (e.g. avoid 3+ significant
## figures on the sd; round generously and say so).

## ============================================================
## Scenario 3: SPATIOTEMPORAL block CV
## Spatial blocking: cluster the 18 city coordinates into K groups
## (k-means on Lat/Lon), hold out one cluster of cities ENTIRELY.
## Temporal blocking: within the retained cities, still respect
## the 2015-2021 / 2022-2024 split.
## This is the strongest test: unseen locations AND unseen future time.
## ============================================================
city_coords <- df %>% distinct(City, Lat, Lon)

K <- 4  # number of spatial clusters; adjust based on how many cities you can afford to hold out
set.seed(20240715)
km <- kmeans(city_coords[, c("Lon", "Lat")], centers = K, nstart = 50)
city_coords$cluster <- km$cluster

cat("\nSpatial clusters (city -> cluster):\n")
print(city_coords %>% arrange(cluster, City))

holdout_cluster <- which.min(table(city_coords$cluster))  # smallest cluster held out
holdout_cities  <- city_coords %>% filter(cluster == holdout_cluster) %>% pull(City)
cat("\nHeld-out spatial cluster:", holdout_cluster,
    "| Held-out cities:", paste(holdout_cities, collapse = ", "), "\n")

df2 <- df %>% left_join(city_coords %>% select(City, cluster), by = "City")

train_st <- df2 %>% filter(!(City %in% holdout_cities), Year <= 2021)
test_st  <- df2 %>% filter(City %in% holdout_cities, Year >= 2022)

results$spatiotemporal <- fit_predict_gwr(train_st, test_st,
                                           "Spatiotemporal block (held-out cities x 2022-24)")

## ---- Optional 4th check: PURE spatial block (all years, held-out cities) ----
## Isolates spatial generalization only, useful to decompose where error comes from
train_spatial_only <- df2 %>% filter(!(City %in% holdout_cities))
test_spatial_only  <- df2 %>% filter(City %in% holdout_cities)

results$spatial_only <- fit_predict_gwr(train_spatial_only, test_spatial_only,
                                         "Spatial block only (held-out cities, all years)")

## ============================================================
## Scenario 5: Leave-One-City-Out (LOCO) spatial block CV
## Runs the spatial holdout 18 times (once per city) instead of once,
## so you get a DISTRIBUTION of RMSE/MAE across independent locations,
## not a single number from one arbitrarily-chosen city.
## ============================================================
all_cities <- sort(unique(df$City))

## Same cluster `cl` created earlier (before Scenario 1) is reused here --
## each city-fold is independent of every other, so this is again
## embarrassingly parallel. df and fit_predict_gwr were already exported
## to the workers; Year/City filtering happens inside the worker itself.
run_loco <- function(spatiotemporal = FALSE) {

  loco_worker <- function(held_city) {
    if (spatiotemporal) {
      train_fold <- df %>% filter(City != held_city, Year <= 2021)
      test_fold  <- df %>% filter(City == held_city, Year >= 2022)
      lbl <- paste0("ST-LOCO city=", held_city)
    } else {
      train_fold <- df %>% filter(City != held_city)
      test_fold  <- df %>% filter(City == held_city)
      lbl <- paste0("Spatial-LOCO city=", held_city)
    }

    ## bw.gwr can occasionally fail to converge on a reduced training set;
    ## don't let one bad fold silently kill the whole run -- flag it instead
    res <- tryCatch(
      fit_predict_gwr(train_fold, test_fold, lbl),
      error = function(e) {
        list(label = lbl, rmse = NA_real_, mae = NA_real_, lag1 = NA_real_,
             n_train = nrow(train_fold), n_test = nrow(test_fold),
             error_msg = conditionMessage(e))
      }
    )
    res$City <- held_city
    res
  }

  ## clusterExport needs `spatiotemporal` (the flag) visible to workers too,
  ## since loco_worker() is defined fresh inside this function each call
  clusterExport(cl, varlist = c("spatiotemporal"), envir = environment())
  fold_results <- parLapply(cl, all_cities, loco_worker)

  failed <- Filter(function(x) is.na(x$rmse), fold_results)
  if (length(failed) > 0) {
    for (fr in failed) {
      message(sprintf("  [FOLD FAILED] city=%s : %s", fr$City,
                       if (!is.null(fr$error_msg)) fr$error_msg else "unknown error"))
    }
  }

  bind_rows(lapply(fold_results, function(x) {
    data.frame(City = x$City, n_train = x$n_train, n_test = x$n_test,
               RMSE = x$rmse, MAE = x$mae, Lag1 = x$lag1)
  }))
}


cat("\n\n==================== SPATIAL LOCO (all years) ====================\n")
loco_spatial <- run_loco(spatiotemporal = FALSE)
print(loco_spatial, row.names = FALSE)

n_failed_spatial <- sum(is.na(loco_spatial$RMSE))
if (n_failed_spatial > 0) {
  warning(sprintf("%d of %d spatial-LOCO folds failed to fit -- excluded from summary stats below.",
                   n_failed_spatial, nrow(loco_spatial)))
}

cat("\n\n==================== SPATIOTEMPORAL LOCO (train<=2021, test>=2022) ====================\n")
loco_st <- run_loco(spatiotemporal = TRUE)
print(loco_st, row.names = FALSE)

n_failed_st <- sum(is.na(loco_st$RMSE))
if (n_failed_st > 0) {
  warning(sprintf("%d of %d spatiotemporal-LOCO folds failed to fit -- excluded from summary stats below.",
                   n_failed_st, nrow(loco_st)))
}

## ---- Fold-level summary stats (mean, sd, median, min, max, n_valid) ----
summarise_loco <- function(loco_df, label) {
  valid <- loco_df %>% filter(!is.na(RMSE))
  data.frame(
    Scenario   = label,
    n_folds    = nrow(loco_df),
    n_valid    = nrow(valid),
    RMSE_mean  = round(mean(valid$RMSE), 4),
    RMSE_sd    = round(sd(valid$RMSE), 4),
    RMSE_median= round(median(valid$RMSE), 4),
    RMSE_min   = round(min(valid$RMSE), 4),
    RMSE_max   = round(max(valid$RMSE), 4),
    MAE_mean   = round(mean(valid$MAE), 4),
    MAE_sd     = round(sd(valid$MAE), 4)
  )
}

loco_summary <- bind_rows(
  summarise_loco(loco_spatial, "Spatial LOCO (18 folds, all years)"),
  summarise_loco(loco_st,      "Spatiotemporal LOCO (18 folds, train<=2021/test>=2022)")
)

cat("\n\n==================== LOCO SUMMARY (mean +/- sd across cities) ====================\n")
print(loco_summary, row.names = FALSE)

write.csv(loco_spatial,  "loco_spatial_by_city.csv",       row.names = FALSE)
write.csv(loco_st,       "loco_spatiotemporal_by_city.csv", row.names = FALSE)
write.csv(loco_summary,  "loco_summary_stats.csv",          row.names = FALSE)
cat("\nSaved: loco_spatial_by_city.csv, loco_spatiotemporal_by_city.csv, loco_summary_stats.csv\n")

## ---- Optional: boxplot of per-city RMSE for the paper's robustness figure ----
## png("loco_rmse_boxplot.png", width = 900, height = 600, res = 120)
## boxplot(loco_spatial$RMSE, loco_st$RMSE,
##         names = c("Spatial LOCO", "Spatiotemporal LOCO"),
##         ylab = "RMSE (per held-out city)",
##         main = "Distribution of RMSE across 18 leave-one-city-out folds")
## dev.off()

## ============================================================
## Summary table for the paper
## ============================================================
summary_tbl <- bind_rows(lapply(results, function(x) {
  data.frame(Scenario = x$label, n_train = x$n_train, n_test = x$n_test,
             RMSE = round(x$rmse, 4),
             RMSE_sd = if (!is.null(x$rmse_sd)) round(x$rmse_sd, 4) else NA_real_,
             MAE = round(x$mae, 4),
             MAE_sd = if (!is.null(x$mae_sd)) round(x$mae_sd, 4) else NA_real_,
             Mean_Lag1_Resid_Autocorr = round(x$lag1, 3))
}))

cat("\n\n==================== SUMMARY ====================\n")
print(summary_tbl, row.names = FALSE)
cat("\nNOTE: 'Random split' RMSE_sd/MAE_sd are the sd across", N_REPS, "repeated 70/30\n")
cat("splits (different seeds) -- see random_split_100reps.csv for the per-rep values.\n")
cat("'Temporal block' is still a single fixed partition (2015-21/2022-24), so it has\n")
cat("no sd -- there is only one way to split by that fixed year boundary.\n")
cat("'Spatiotemporal block' and 'Spatial block only' rows are from ONE arbitrarily-\n")
cat("chosen held-out city (the smallest k-means cluster) and have no sd either; for\n")
cat("the paper, prefer the LOCO mean +/- sd numbers in loco_summary_stats.csv, which\n")
cat("are based on 18 independent folds and are what a reviewer will trust more.\n")

write.csv(summary_tbl, "blockcv_comparison_results_single_fold.csv", row.names = FALSE)
cat("\nSaved: blockcv_comparison_results_single_fold.csv\n")

## ============================================================
## Clean up the parallel cluster explicitly at the end too (in addition
## to the on.exit() safety net registered earlier)
## ============================================================
stopCluster(cl)
cat("\nParallel cluster stopped. Done.\n")
