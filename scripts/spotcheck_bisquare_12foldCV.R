## ============================================================
## *** ARCHIVAL NOTE (added after unifying the random-split ratio to
## 75/25 across the paper; see rerun_sr_random_split_75_25.R) ***
##
## The "Random split" section below (Spot-check 1) uses a 70/30 split
## (floor(0.7 * n_total)) and is SUPERSEDED for the numbers reported in
## the final paper. The bi-square Random-split baseline actually
## reported in the paper (RMSE = 0.6785 +/- 0.0267) was produced by
## rerun_sr_random_split_75_25.R using a 75/25 split instead; use that
## script's output (spotcheck_random_bisquare.csv) as the authoritative
## source for that number.
##
## The "Spatial LOCO" section below (Spot-check 2) does NOT depend on
## the random-split ratio and remains fully valid and unchanged; its
## output (spotcheck_spatial_loco_bisquare.csv) is the authoritative
## source for the bi-square Spatial LOCO numbers used throughout the
## paper.
## ============================================================

## ============================================================
## Spot-check: bi-square kernel + manual 12-fold CV bandwidth selection,
## matching the SR configuration used in Table 2 (paper's main results),
## instead of the Gaussian + AICc configuration used in the robustness
## analysis (Section 3.4 / gwr_blockcv_scenarios.R).
##
## PURPOSE: verify that the robustness conclusions (spatial LOCO >>
## random split; temporal leakage negligible) are not artifacts of the
## AICc/Gaussian shortcut used for computational speed. This script does
## NOT redo the full 141-fold analysis -- it only checks two scenarios:
## (1) Random split, 20 repetitions (2) Spatial LOCO, all 18 folds.
##
## RUNTIME WARNING: this is much slower per fold than the AICc version,
## because bandwidth selection itself now requires a nested 12-fold CV
## over a grid of candidate bandwidths. Expect this to take substantially
## longer than the original script, even in parallel. Consider running
## it unattended (overnight) or on an HPC/Slurm allocation if available.
## ============================================================

set.seed(20240715)

required_pkgs <- c("readxl", "dplyr", "GWmodel", "sp", "parallel")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

PROJECT_ROOT <- Sys.getenv("IRAQ_WEATHER_ROOT", unset = getwd())  # set by run_all.R; see README
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "03_spotcheck_bisquare_12foldCV")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)

library(readxl); library(dplyr); library(GWmodel); library(sp); library(parallel)

rm(list = setdiff(ls(), c("PROJECT_ROOT", "OUTPUT_DIR")))
set.seed(20240715)

## ---- 1. Load & clean ----
path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
raw <- read_excel(path)  ## read the data ("Dalya5x.xlsx")


#raw <- read_excel("Dalya5x.xlsx")
df <- raw %>%
  rename(City = City, Humidity = Humidity, Temp = `Temp.`,
         TOA = `Top of Atmosphere`, Wind = `Wind Speed`,
         Precip = Precipitation, Soil = `Soil Moisture`,
         Lat = Latitude, Lon = Longitude) %>%
  arrange(City) %>%
  group_by(City) %>%
  mutate(month_idx = row_number(),
         Year = 2015 + (month_idx - 1) %/% 12,
         Month = ((month_idx - 1) %% 12) + 1) %>%
  ungroup()

model_formula <- Humidity ~ Temp + TOA + Wind + Precip + Soil

## ============================================================
## Manual k-fold CV bandwidth selection (bi-square kernel, adaptive)
## Mirrors what the paper's Methods section describes for the SR model
## in Table 2, since GWmodel's built-in bw.gwr() only offers AICc or
## leave-one-out CV, not k-fold CV.
## ============================================================
select_bw_kfoldCV <- function(train_df, k = 12, bw_grid = NULL, seed = 1) {
  n <- nrow(train_df)

  if (is.null(bw_grid)) {
    ## FIX: the data has only 18 distinct spatial locations, each repeated
    ## many times over (monthly records). Coordinates are literally
    ## identical for all rows of the same city, so the nearest neighbors
    ## of any point are always same-city points first. If the adaptive
    ## bandwidth (neighbor count) is smaller than a city's own sample
    ## size, the "spatial" neighborhood collapses onto a single city's
    ## time series -- a degenerate local fit that produces zero-weight
    ## (NaN) predictions for test points in other cities. We floor the
    ## candidate grid safely above the largest per-city count actually
    ## present in this training fold, so every candidate bandwidth is
    ## guaranteed to span at least two cities.
    max_city_n <- train_df %>% count(City) %>% pull(n) %>% max()
    floor_bw <- max_city_n + 10  # a small safety margin above one full city
    bw_grid <- unique(round(seq(floor_bw, min(n - 1, 500), length.out = 15)))
  }

  set.seed(seed)
  folds <- sample(rep(1:k, length.out = n))

  cv_mse_for_bw <- function(bw) {
    fold_errors <- numeric(k)
    for (fold in seq_len(k)) {
      tr <- train_df[folds != fold, ]
      te <- train_df[folds == fold, ]

      tr_sp <- tr; coordinates(tr_sp) <- ~ Lon + Lat
      te_sp <- te; coordinates(te_sp) <- ~ Lon + Lat

      pred <- tryCatch(
        gwr.predict(model_formula, data = tr_sp, predictdata = te_sp,
                    bw = bw, kernel = "bisquare", adaptive = TRUE)$SDF$prediction,
        error = function(e) rep(NA_real_, nrow(te))
      )

      sq_err <- (te$Humidity - pred)^2
      ## FIX: do NOT silently na.rm non-finite errors -- a bandwidth that
      ## produces even one undefined (NaN/Inf) prediction is unstable and
      ## must be penalized, not scored on the "easy" subset that happened
      ## to work. This is what let bw=64 look artificially good before.
      if (any(!is.finite(sq_err))) {
        fold_errors[fold] <- Inf
      } else {
        fold_errors[fold] <- mean(sq_err)
      }
    }
    mean(fold_errors)  # if any fold was Inf, this bw's score is Inf -- never selected
  }

  cv_scores <- sapply(bw_grid, cv_mse_for_bw)

  if (all(!is.finite(cv_scores))) {
    stop("All candidate bandwidths produced unstable (non-finite) predictions ",
         "in at least one CV fold -- widen or shift bw_grid upward and retry.")
  }

  best_bw <- bw_grid[which.min(cv_scores)]
  list(best_bw = best_bw, bw_grid = bw_grid, cv_scores = cv_scores)
}

## ============================================================
## Fit/predict wrapper using bi-square kernel + 12-fold CV bandwidth
## (drop-in replacement for the Gaussian/AICc fit_predict_gwr used
## elsewhere, but MUCH slower per call -- see runtime warning above)
## ============================================================
fit_predict_gwr_bisquare <- function(train_df, test_df, label) {
  bw_sel <- select_bw_kfoldCV(train_df, k = 12)
  bw_best <- bw_sel$best_bw

  train_sp <- train_df; coordinates(train_sp) <- ~ Lon + Lat
  test_sp  <- test_df;  coordinates(test_sp)  <- ~ Lon + Lat

  pred <- gwr.predict(model_formula, data = train_sp, predictdata = test_sp,
                       bw = bw_best, kernel = "bisquare", adaptive = TRUE)

  yhat <- pred$SDF$prediction
  yobs <- test_df$Humidity
  resid <- yobs - yhat
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  mae  <- mean(abs(resid), na.rm = TRUE)

  cat(sprintf("[%s] bw(bisquare,12-fold CV) = %d | RMSE = %.4f | MAE = %.4f\n",
              label, bw_best, rmse, mae))

  list(label = label, rmse = rmse, mae = mae, bw = bw_best,
       n_train = nrow(train_df), n_test = nrow(test_df))
}

## ============================================================
## Spot-check 1: Random split, 20 repetitions (not 100 -- this
## configuration is far more expensive per fold, so 20 reps is a
## reasonable compromise; note this explicitly if you use it in the paper)
##
## CHECKPOINTING: results are appended to disk after EVERY repetition,
## not just at the end. If the run is interrupted for any reason, you
## can resume from the CSV instead of losing everything computed so far.
## ============================================================
N_REPS_SPOTCHECK <- 20
n_total <- nrow(df)

checkpoint_file_random <- "spotcheck_random_bisquare.csv"

## Resume support: skip reps already completed in a previous run
already_done_random <- integer(0)
if (file.exists(checkpoint_file_random)) {
  prev <- read.csv(checkpoint_file_random)
  ## FIX: only treat a rep as "done" if it actually succeeded (finite RMSE).
  ## Reps that produced NaN under the old bw_grid (e.g. bw=64) will be
  ## retried automatically now that the bandwidth floor has been fixed.
  if ("rep_id" %in% names(prev)) {
    already_done_random <- prev$rep_id[is.finite(prev$RMSE)]
  }
  n_failed_prev <- sum(!is.finite(prev$RMSE))
  cat(sprintf("Resuming: %d reps already succeeded, %d will be retried with the fix\n",
              length(already_done_random), n_failed_prev))
} else {
  ## write header once. NOTE: column renamed rep -> rep_id, because `rep`
  ## is also the name of a base R function (used for repeating vectors);
  ## a column literally named `rep` previously caused dplyr::arrange(rep, ...)
  ## to resolve to the FUNCTION rep() instead of the column, crashing with
  ## "must be a vector, not a primitive function". Same class of bug as the
  ## earlier `n` vs. dplyr::n() collision -- avoid single-letter / common
  ## function names for data.frame columns whenever dplyr is loaded.
  write.csv(data.frame(rep_id = integer(0), RMSE = numeric(0), MAE = numeric(0), bw = integer(0)),
            checkpoint_file_random, row.names = FALSE)
}

for (r in seq_len(N_REPS_SPOTCHECK)) {
  if (r %in% already_done_random) next  # skip work already done

  set.seed(20240715 + r)
  rand_idx <- sample.int(n_total, size = floor(0.7 * n_total))
  train_rand <- df[rand_idx, ]
  test_rand  <- df[-rand_idx, ]

  res <- tryCatch(
    fit_predict_gwr_bisquare(train_rand, test_rand, sprintf("Random(bisquare) rep=%d", r)),
    error = function(e) {
      message(sprintf("  [FAILED] rep=%d : %s", r, conditionMessage(e)))
      list(label = sprintf("rep=%d", r), rmse = NA_real_, mae = NA_real_, bw = NA_integer_)
    }
  )

  ## Append this single result immediately -- this line is the fix.
  write.table(data.frame(rep_id = r, RMSE = res$rmse, MAE = res$mae, bw = res$bw),
              checkpoint_file_random, sep = ",", row.names = FALSE, col.names = FALSE,
              append = TRUE)
  cat(sprintf("  [checkpoint saved] rep=%d/%d\n", r, N_REPS_SPOTCHECK))
}

random_bisquare_df <- read.csv(checkpoint_file_random)

## SANITY FILTER: drop any row with an implausible RMSE. For this humidity
## variable (range ~2.6-9.6 g/kg), a valid RMSE should be well under 5; any
## row far outside that (e.g. from an earlier corrupted/merged checkpoint
## file) is discarded rather than silently averaged in.
n_before <- nrow(random_bisquare_df)
random_bisquare_df <- random_bisquare_df[is.finite(random_bisquare_df$RMSE) &
                                          random_bisquare_df$RMSE > 0 &
                                          random_bisquare_df$RMSE < 5, ]
if (nrow(random_bisquare_df) < n_before) {
  warning(sprintf("Dropped %d implausible/corrupted row(s) from %s -- inspect the raw file if this number is large.",
                   n_before - nrow(random_bisquare_df), checkpoint_file_random))
}

## Deduplicate with base R only (no dplyr::arrange/distinct -- avoids any
## further risk of column-name/function-name masking): for each rep_id,
## keep the first row encountered (they are all valid at this point, since
## the sanity filter above already removed non-finite/implausible rows).
random_bisquare_df <- random_bisquare_df[!duplicated(random_bisquare_df$rep_id), ]
write.csv(random_bisquare_df, checkpoint_file_random, row.names = FALSE)

cat(sprintf("\nRandom split (bi-square, %d-fold CV bw), %d reps: RMSE = %.4f +/- %.4f\n",
            12, nrow(random_bisquare_df),
            mean(random_bisquare_df$RMSE, na.rm = TRUE), sd(random_bisquare_df$RMSE, na.rm = TRUE)))

## ============================================================
## Spot-check 2: Spatial LOCO, all 18 folds
## Same checkpointing pattern: append after every city, resume-safe.
## ============================================================
all_cities <- sort(unique(df$City))
checkpoint_file_spatial <- "spotcheck_spatial_loco_bisquare.csv"

already_done_spatial <- integer(0)
if (file.exists(checkpoint_file_spatial)) {
  prev <- read.csv(checkpoint_file_spatial)
  if ("City" %in% names(prev)) {
    already_done_spatial <- prev$City[is.finite(prev$RMSE)]
  }
  n_failed_prev <- sum(!is.finite(prev$RMSE))
  cat(sprintf("Resuming: %d cities already succeeded, %d will be retried with the fix\n",
              length(already_done_spatial), n_failed_prev))
} else {
  write.csv(data.frame(City = integer(0), RMSE = numeric(0), MAE = numeric(0), bw = integer(0)),
            checkpoint_file_spatial, row.names = FALSE)
}

for (i in seq_along(all_cities)) {
  held_city <- all_cities[i]
  if (held_city %in% already_done_spatial) next

  train_fold <- df %>% filter(City != held_city)
  test_fold  <- df %>% filter(City == held_city)

  res <- tryCatch(
    fit_predict_gwr_bisquare(train_fold, test_fold, paste0("Spatial-LOCO(bisquare) city=", held_city)),
    error = function(e) {
      message(sprintf("  [FAILED] city=%s : %s", held_city, conditionMessage(e)))
      list(label = paste0("city=", held_city), rmse = NA_real_, mae = NA_real_, bw = NA_integer_)
    }
  )

  write.table(data.frame(City = held_city, RMSE = res$rmse, MAE = res$mae, bw = res$bw),
              checkpoint_file_spatial, sep = ",", row.names = FALSE, col.names = FALSE,
              append = TRUE)
  cat(sprintf("  [checkpoint saved] city=%s (%d/%d)\n", held_city, i, length(all_cities)))
}

spatial_bisquare_df <- read.csv(checkpoint_file_spatial)

n_before_sp <- nrow(spatial_bisquare_df)
spatial_bisquare_df <- spatial_bisquare_df[is.finite(spatial_bisquare_df$RMSE) &
                                            spatial_bisquare_df$RMSE > 0 &
                                            spatial_bisquare_df$RMSE < 5, ]
if (nrow(spatial_bisquare_df) < n_before_sp) {
  warning(sprintf("Dropped %d implausible/corrupted row(s) from %s.",
                   n_before_sp - nrow(spatial_bisquare_df), checkpoint_file_spatial))
}

spatial_bisquare_df <- spatial_bisquare_df[!duplicated(spatial_bisquare_df$City), ]
write.csv(spatial_bisquare_df, checkpoint_file_spatial, row.names = FALSE)

cat(sprintf("\nSpatial LOCO (bi-square, %d-fold CV bw), %d folds: RMSE = %.4f +/- %.4f\n",
            12, nrow(spatial_bisquare_df),
            mean(spatial_bisquare_df$RMSE, na.rm = TRUE), sd(spatial_bisquare_df$RMSE, na.rm = TRUE)))

## ============================================================
## Quick significance check: does Spatial LOCO (bi-square) still
## differ significantly from the Random split (bi-square) null?
## ============================================================
wt <- wilcox.test(spatial_bisquare_df$RMSE, random_bisquare_df$RMSE, alternative = "greater")
print(wt)

cat(sprintf("\nMann-Whitney (Spatial LOCO > Random, bi-square config): p = %.5f\n", wt$p.value))
cat("If this p-value is still small (e.g. < 0.05), the paper's main robustness\n")
cat("conclusion (spatial leakage matters, temporal leakage doesn't) is NOT an\n")
cat("artifact of the Gaussian/AICc shortcut, and you can state this explicitly\n")
cat("in the footnote instead of leaving it as an open caveat.\n")

