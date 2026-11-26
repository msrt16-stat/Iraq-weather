
PROJECT_ROOT <- Sys.getenv("IRAQ_WEATHER_ROOT", unset = getwd())  # set by run_all.R; see README
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "04_rerun_sr_random_split_75_25")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)
getwd()



## ============================================================
## Re-run ONLY the SR random-split baselines with rTrain=0.75 (75/25),
## to unify the split ratio with Table 2 across the entire robustness
## analysis in Section "Robustness to Cross-Validation Design".
##
## WHY ONLY THESE TWO PIECES: the LOCO, spatiotemporal LOCO, and
## rolling-origin folds do not depend on a random train/test percentage
## at all (they split by governorate or by year), so they are UNCHANGED
## and do not need to be re-run. Only the two "repeated random split"
## baselines (Gaussian/AICc config and bi-square/12-fold-CV config) used
## 70/30 and need to be redone at 75/25.
##
## OUTPUT: this script overwrites two files with the SAME names used
## throughout the rest of the analysis, so no other script or number
## needs to change:
##   random_split_100reps.csv           (Gaussian/AICc config, 100 reps)
##   spotcheck_random_bisquare.csv      (bi-square/12-fold-CV config, 20 reps)
## It then recomputes the Mann-Whitney comparisons against the EXISTING,
## unchanged LOCO/rolling-origin CSVs and prints updated test statistics
## and effect sizes for you to paste into the paper.
##
## RUNTIME: substantially cheaper than the original full analysis, since
## only the random-split reps are redone, not any LOCO fold.
## ============================================================

library(openxlsx)
library(spgwr)
library(sf)
library(sp)
library(GWmodel)

set.seed(20240715)
rTrain <- 0.75   # <-- the only real change: was 0.7 in the original SR robustness scripts

path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
raw <- read.xlsx(path)
MyData <- raw
colnames(MyData) <- c("City","Humidity","Temp","TOA","Wind","Precip","Soil","Lat","Lon")
MyData <- MyData[order(MyData$City), ]
MyData$month_idx <- ave(MyData$City, MyData$City, FUN = seq_along)
MyData$Year <- 2015 + (MyData$month_idx - 1) %/% 12

model_formula <- Humidity ~ Temp + TOA + Wind + Precip + Soil

## ============================================================
## Config A: Gaussian kernel, AICc bandwidth (bw.gwr default) --
## matches the original random_split_100reps.csv analysis
## ============================================================
fit_predict_gwr_gaussian <- function(train_df, test_df) {
  train_sp <- train_df; coordinates(train_sp) <- ~ Lon + Lat
  test_sp  <- test_df;  coordinates(test_sp)  <- ~ Lon + Lat

  bw_train <- bw.gwr(model_formula, data = train_sp, approach = "AICc",
                      kernel = "gaussian", adaptive = TRUE)
  pred <- gwr.predict(model_formula, data = train_sp, predictdata = test_sp,
                       bw = bw_train, kernel = "gaussian", adaptive = TRUE)

  resid <- test_df$Humidity - pred$SDF$prediction
  list(rmse = sqrt(mean(resid^2, na.rm = TRUE)), mae = mean(abs(resid), na.rm = TRUE))
}

N_REPS_GAUSSIAN <- 100
checkpoint_gaussian <- "random_split_100reps.csv"

already_done_g <- integer(0)
if (file.exists(checkpoint_gaussian)) {
  prev <- read.csv(checkpoint_gaussian)
  if ("rep_id" %in% names(prev)) already_done_g <- prev$rep_id[is.finite(prev$RMSE)]
  cat(sprintf("Resuming Gaussian config: %d/%d reps already done\n", length(already_done_g), N_REPS_GAUSSIAN))
} else {
  write.csv(data.frame(rep_id = integer(0), RMSE = numeric(0), MAE = numeric(0)),
            checkpoint_gaussian, row.names = FALSE)
}

n_total <- nrow(MyData)
for (r in seq_len(N_REPS_GAUSSIAN)) {
  if (r %in% already_done_g) next
  set.seed(20240715 + r)
  idx <- sample(seq_len(n_total), size = round(rTrain * n_total))
  train_r <- MyData[idx, ]; test_r <- MyData[-idx, ]

  res <- tryCatch(fit_predict_gwr_gaussian(train_r, test_r),
                   error = function(e) { message(sprintf("  [FAILED] rep=%d: %s", r, conditionMessage(e)))
                                          list(rmse = NA_real_, mae = NA_real_) })
  cat(sprintf("[Gaussian, 75/25] rep=%d: RMSE=%.4f MAE=%.4f\n", r, res$rmse, res$mae))
  write.table(data.frame(rep_id = r, RMSE = res$rmse, MAE = res$mae),
              checkpoint_gaussian, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

gaussian_results <- read.csv(checkpoint_gaussian)
gaussian_results <- gaussian_results[is.finite(gaussian_results$RMSE), ]
gaussian_results <- gaussian_results[!duplicated(gaussian_results$rep_id), ]
write.csv(gaussian_results, checkpoint_gaussian, row.names = FALSE)
cat(sprintf("\nGaussian/AICc random split (75/25, n=%d): RMSE = %.4f +/- %.4f\n\n",
            nrow(gaussian_results), mean(gaussian_results$RMSE), sd(gaussian_results$RMSE)))

## ============================================================
## Config B: bi-square kernel, 12-fold CV bandwidth selection --
## matches the original spotcheck_random_bisquare.csv analysis.
## Includes the governorate-size bandwidth floor fix established earlier
## (avoids the degenerate small-bandwidth / NaN issue).
## ============================================================
select_bw_kfoldCV <- function(train_df, k = 12, bw_grid = NULL, seed = 1) {
  n <- nrow(train_df)
  if (is.null(bw_grid)) {
    max_city_n <- max(table(train_df$City))
    floor_bw <- max_city_n + 10
    bw_grid <- unique(round(seq(floor_bw, min(n - 1, 500), length.out = 15)))
  }
  set.seed(seed)
  folds <- sample(rep(1:k, length.out = n))

  cv_mse_for_bw <- function(bw) {
    fold_errors <- numeric(k)
    for (fold in seq_len(k)) {
      tr <- train_df[folds != fold, ]; te <- train_df[folds == fold, ]
      tr_sp <- tr; coordinates(tr_sp) <- ~ Lon + Lat
      te_sp <- te; coordinates(te_sp) <- ~ Lon + Lat
      pred <- tryCatch(gwr.predict(model_formula, data = tr_sp, predictdata = te_sp,
                                    bw = bw, kernel = "bisquare", adaptive = TRUE)$SDF$prediction,
                        error = function(e) rep(NA_real_, nrow(te)))
      sq_err <- (te$Humidity - pred)^2
      fold_errors[fold] <- if (any(!is.finite(sq_err))) Inf else mean(sq_err)
    }
    mean(fold_errors)
  }
  cv_scores <- sapply(bw_grid, cv_mse_for_bw)
  if (all(!is.finite(cv_scores))) stop("All candidate bandwidths unstable -- widen bw_grid.")
  list(best_bw = bw_grid[which.min(cv_scores)])
}

fit_predict_gwr_bisquare <- function(train_df, test_df) {
  bw_sel <- select_bw_kfoldCV(train_df, k = 12)
  train_sp <- train_df; coordinates(train_sp) <- ~ Lon + Lat
  test_sp  <- test_df;  coordinates(test_sp)  <- ~ Lon + Lat
  pred <- gwr.predict(model_formula, data = train_sp, predictdata = test_sp,
                       bw = bw_sel$best_bw, kernel = "bisquare", adaptive = TRUE)
  resid <- test_df$Humidity - pred$SDF$prediction
  list(rmse = sqrt(mean(resid^2, na.rm = TRUE)), mae = mean(abs(resid), na.rm = TRUE),
       bw = bw_sel$best_bw)
}

N_REPS_BISQUARE <- 20
checkpoint_bisquare <- "spotcheck_random_bisquare.csv"

already_done_b <- integer(0)
if (file.exists(checkpoint_bisquare)) {
  prev <- read.csv(checkpoint_bisquare)
  if ("rep_id" %in% names(prev)) already_done_b <- prev$rep_id[is.finite(prev$RMSE)]
  cat(sprintf("Resuming bi-square config: %d/%d reps already done\n", length(already_done_b), N_REPS_BISQUARE))
} else {
  write.csv(data.frame(rep_id = integer(0), RMSE = numeric(0), MAE = numeric(0), bw = integer(0)),
            checkpoint_bisquare, row.names = FALSE)
}

for (r in seq_len(N_REPS_BISQUARE)) {
  if (r %in% already_done_b) next
  set.seed(20240715 + r)
  idx <- sample(seq_len(n_total), size = round(rTrain * n_total))
  train_r <- MyData[idx, ]; test_r <- MyData[-idx, ]

  res <- tryCatch(fit_predict_gwr_bisquare(train_r, test_r),
                   error = function(e) { message(sprintf("  [FAILED] rep=%d: %s", r, conditionMessage(e)))
                                          list(rmse = NA_real_, mae = NA_real_, bw = NA_integer_) })
  cat(sprintf("[bi-square, 75/25] rep=%d: RMSE=%.4f MAE=%.4f bw=%s\n", r, res$rmse, res$mae, res$bw))
  write.table(data.frame(rep_id = r, RMSE = res$rmse, MAE = res$mae, bw = res$bw),
              checkpoint_bisquare, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

bisquare_results <- read.csv(checkpoint_bisquare)
bisquare_results <- bisquare_results[is.finite(bisquare_results$RMSE), ]
bisquare_results <- bisquare_results[!duplicated(bisquare_results$rep_id), ]
write.csv(bisquare_results, checkpoint_bisquare, row.names = FALSE)
cat(sprintf("\nBi-square/12-fold-CV random split (75/25, n=%d): RMSE = %.4f +/- %.4f\n\n",
            nrow(bisquare_results), mean(bisquare_results$RMSE), sd(bisquare_results$RMSE)))

## ============================================================
## Recompute all Mann-Whitney comparisons that use these baselines,
## against the EXISTING (unchanged) LOCO / rolling-origin CSVs.
## Adjust these file paths if your LOCO CSVs are named differently.
## ============================================================
mannwhitney_report <- function(high, low, label) {
  wt <- wilcox.test(high, low, alternative = "greater")
  n1 <- length(high); n2 <- length(low)
  cles <- wt$statistic / (n1 * n2)
  r <- 2 * cles - 1
  cat(sprintf("%s: U=%.1f p=%.4g r=%.3f (n1=%d, n2=%d)\n", label, wt$statistic, wt$p.value, r, n1, n2))
}

loco_sp_path      <- file.path(PROJECT_ROOT, "outputs", "02_gwr_blockcv_scenarios", "loco_spatial_by_city.csv")
loco_st_path      <- file.path(PROJECT_ROOT, "outputs", "02_gwr_blockcv_scenarios", "loco_spatiotemporal_by_city.csv")
rolling_path      <- file.path(PROJECT_ROOT, "outputs", "02_gwr_blockcv_scenarios", "rolling_temporal_by_year.csv")
bisq_loco_path    <- file.path(PROJECT_ROOT, "outputs", "03_spotcheck_bisquare_12foldCV", "spotcheck_spatial_loco_bisquare.csv")

if (file.exists(loco_sp_path)) {
  loco_sp <- read.csv(loco_sp_path)$RMSE
  mannwhitney_report(loco_sp, gaussian_results$RMSE, "Spatial LOCO (Gaussian) vs Random (75/25)")
} else {
  cat(sprintf("NOTE: %s not found -- run 02_gwr_blockcv_scenarios first.\n", loco_sp_path))
}
if (file.exists(loco_st_path)) {
  loco_st <- read.csv(loco_st_path)$RMSE
  mannwhitney_report(loco_st, gaussian_results$RMSE, "Spatiotemporal LOCO (Gaussian) vs Random (75/25)")
} else {
  cat(sprintf("NOTE: %s not found -- run 02_gwr_blockcv_scenarios first.\n", loco_st_path))
}
if (file.exists(rolling_path)) {
  rolling <- read.csv(rolling_path)$RMSE
  mannwhitney_report(rolling, gaussian_results$RMSE, "Rolling-origin vs Random (75/25)")
} else {
  cat(sprintf("NOTE: %s not found -- run 02_gwr_blockcv_scenarios first.\n", rolling_path))
}
if (file.exists(bisq_loco_path)) {
  bisq_loco <- read.csv(bisq_loco_path)
  mannwhitney_report(bisq_loco$RMSE, bisquare_results$RMSE, "Bi-square Spatial LOCO vs Random (75/25)")
  excl <- bisq_loco[bisq_loco$City != 5, ]$RMSE  # Kirkuk = City ID 5
  mannwhitney_report(excl, bisquare_results$RMSE, "Bi-square Spatial LOCO excl-Kirkuk vs Random (75/25)")
} else {
  cat(sprintf("NOTE: %s not found -- run 03_spotcheck_bisquare_12foldCV first.\n", bisq_loco_path))
}

cat("\nDone. random_split_100reps.csv and spotcheck_random_bisquare.csv now reflect 75/25.\n")

