## ============================================================
## Spatial LOCO validation for the GLOBAL models (plain RF, plain SVR --
## no spatial-lag predictors), to directly test the paper's core claim:
## does incorporating spatial correlation structure (SRF vs RF, SSVR vs
## SVR) improve accuracy and robustness, INCLUDING under genuine spatial
## extrapolation to an unseen city -- not just under random-split
## evaluation where leakage/interpolation could inflate the comparison.
##
## This produces the missing piece needed for a fair, direct comparison:
##   RF-LOCO   vs SRF-LOCO   (already computed in srf_spatial_loco.csv)
##   SVR-LOCO  vs SSVR-LOCO  (already computed in ssvr_spatial_loco.csv)
##
## RUNTIME: cheaper than SRF/SSVR since there is no spatial-lag feature
## to build and the formula has 5 predictors instead of 11, but the RF
## grid search (12-fold CV x 9 ntree/mtry combos) is still the slower
## half of this script.
## ============================================================

PROJECT_ROOT <- "C:/Weather Project"
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "07_global_models_spatial_loco")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)

library(openxlsx)
library(e1071)
library(randomForest)

Seed <- 123

path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
MyData <- read.xlsx(path)
colnames(MyData) <- c("ID","Y","X1","X2","X3","X4","X5","Latitude","Longitude")

Formula <- Y ~ X1 + X2 + X3 + X4 + X5   # no WY/WX -- these are the GLOBAL models

## ============================================================
## RF: full 12-fold CV grid search + fit + predict, for one train/test
## partition (used identically for LOCO folds and random-split reps)
## ============================================================
fit_predict_rf <- function(TrainData, TestData, label,
                            ntree.grid = c(100, 500, 1000), mtry.grid = c(2, 3, 4)) {
  set.seed(Seed)
  K <- 12
  n_k <- nrow(TrainData)
  fold.id <- sample(rep(1:K, length.out = n_k))

  grid.results <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())
  for (nt in ntree.grid) {
    for (mt in mtry.grid) {
      rmse.fold <- numeric(K)
      for (kk in 1:K) {
        Train.CV <- TrainData[fold.id != kk, ]
        Valid.CV <- TrainData[fold.id == kk, ]
        set.seed(Seed)
        rf.cv <- randomForest(Formula, data = Train.CV, ntree = nt, mtry = mt)
        pred.cv <- predict(rf.cv, newdata = Valid.CV)
        rmse.fold[kk] <- sqrt(mean((Valid.CV$Y - pred.cv)^2))
      }
      grid.results <- rbind(grid.results, data.frame(ntree = nt, mtry = mt, RMSE = mean(rmse.fold)))
    }
  }
  best <- grid.results[which.min(grid.results$RMSE), ]

  model <- randomForest(Formula, data = TrainData, ntree = best$ntree, mtry = best$mtry)
  pred  <- predict(model, newdata = TestData)

  rmse <- sqrt(mean((TestData$Y - pred)^2))
  mae  <- mean(abs(TestData$Y - pred))
  r2   <- 1 - sum((TestData$Y - pred)^2) / sum((TestData$Y - mean(TestData$Y))^2)

  cat(sprintf("[%s] R2=%.4f RMSE=%.4f MAE=%.4f (ntree=%d, mtry=%d)\n",
              label, r2, rmse, mae, best$ntree, best$mtry))
  list(R2 = r2, RMSE = rmse, MAE = mae)
}

## ============================================================
## SVR: 12-fold CV tune() + fit + predict
## ============================================================
fit_predict_svr <- function(TrainData, TestData, label) {
  set.seed(Seed)
  tune.svr <- tune(svm, Formula, data = TrainData,
                    type = "eps-regression", kernel = "radial",
                    ranges = list(cost = c(0.1,1,10), gamma = c(0.01,0.1,1)),
                    tunecontrol = tune.control(cross = 12))
  model <- tune.svr$best.model
  pred  <- predict(model, newdata = TestData)

  rmse <- sqrt(mean((TestData$Y - pred)^2))
  mae  <- mean(abs(TestData$Y - pred))
  r2   <- 1 - sum((TestData$Y - pred)^2) / sum((TestData$Y - mean(TestData$Y))^2)

  cat(sprintf("[%s] R2=%.4f RMSE=%.4f MAE=%.4f (cost=%.3g, gamma=%.3g)\n",
              label, r2, rmse, mae, tune.svr$best.parameters$cost, tune.svr$best.parameters$gamma))
  list(R2 = r2, RMSE = rmse, MAE = mae)
}

## ============================================================
## Generic runner: LOCO (18 folds) + repeated random-split baseline,
## with checkpointing, for a given fit_predict_* function and file prefix
## ============================================================
run_loco_and_baseline <- function(fit_predict_fn, prefix, n_reps_random = 20) {

  all_cities <- sort(unique(MyData$ID))
  checkpoint_loco <- paste0(prefix, "_spatial_loco.csv")

  already_done <- integer(0)
  if (file.exists(checkpoint_loco)) {
    prev <- read.csv(checkpoint_loco)
    already_done <- prev$City[is.finite(prev$RMSE)]
    cat(sprintf("[%s] Resuming LOCO: %d cities already done\n", prefix, length(already_done)))
  } else {
    write.csv(data.frame(City=integer(0), R2=numeric(0), RMSE=numeric(0), MAE=numeric(0)),
              checkpoint_loco, row.names = FALSE)
  }

  for (city_id in all_cities) {
    if (city_id %in% already_done) next
    TrainData_fold <- MyData[MyData$ID != city_id, ]
    TestData_fold  <- MyData[MyData$ID == city_id, ]

    res <- tryCatch(
      fit_predict_fn(TrainData_fold, TestData_fold, sprintf("%s-LOCO city=%d", prefix, city_id)),
      error = function(e) {
        message(sprintf("  [FAILED] city=%d : %s", city_id, conditionMessage(e)))
        list(R2 = NA, RMSE = NA, MAE = NA)
      }
    )
    write.table(data.frame(City = city_id, R2 = res$R2, RMSE = res$RMSE, MAE = res$MAE),
                checkpoint_loco, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }

  loco_results <- read.csv(checkpoint_loco)
  loco_results <- loco_results[!duplicated(loco_results$City), ]
  write.csv(loco_results, checkpoint_loco, row.names = FALSE)

  ## ---- Repeated random-split baseline ----
  checkpoint_random <- paste0(prefix, "_random_split_baseline.csv")
  already_done_r <- integer(0)
  if (file.exists(checkpoint_random)) {
    prev_r <- read.csv(checkpoint_random)
    already_done_r <- prev_r$rep_id[is.finite(prev_r$RMSE)]
    cat(sprintf("[%s] Resuming random baseline: %d reps already done\n", prefix, length(already_done_r)))
  } else {
    write.csv(data.frame(rep_id=integer(0), R2=numeric(0), RMSE=numeric(0), MAE=numeric(0)),
              checkpoint_random, row.names = FALSE)
  }

  n_total <- nrow(MyData)
  for (r in seq_len(n_reps_random)) {
    if (r %in% already_done_r) next
    set.seed(20240715 + r)
    idx <- sample(1:n_total, size = round(0.75 * n_total))
    TrainData_r <- MyData[idx, ]
    TestData_r  <- MyData[-idx, ]

    res_r <- tryCatch(
      fit_predict_fn(TrainData_r, TestData_r, sprintf("%s-Random rep=%d", prefix, r)),
      error = function(e) {
        message(sprintf("  [FAILED] rep=%d : %s", r, conditionMessage(e)))
        list(R2 = NA, RMSE = NA, MAE = NA)
      }
    )
    write.table(data.frame(rep_id = r, R2 = res_r$R2, RMSE = res_r$RMSE, MAE = res_r$MAE),
                checkpoint_random, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }

  random_results <- read.csv(checkpoint_random)
  random_results <- random_results[!duplicated(random_results$rep_id), ]
  write.csv(random_results, checkpoint_random, row.names = FALSE)

  wt <- wilcox.test(loco_results$RMSE, random_results$RMSE, alternative = "greater")
  delta_pct <- (mean(loco_results$RMSE, na.rm=TRUE) - mean(random_results$RMSE, na.rm=TRUE)) /
               mean(random_results$RMSE, na.rm=TRUE) * 100

  cat(sprintf("\n==================== %s SUMMARY ====================\n", toupper(prefix)))
  cat(sprintf("Random split : RMSE = %.4f +/- %.4f (n=%d)\n",
              mean(random_results$RMSE, na.rm=TRUE), sd(random_results$RMSE, na.rm=TRUE), nrow(random_results)))
  cat(sprintf("Spatial LOCO : RMSE = %.4f +/- %.4f (n=%d)\n",
              mean(loco_results$RMSE, na.rm=TRUE), sd(loco_results$RMSE, na.rm=TRUE), nrow(loco_results)))
  cat(sprintf("Relative increase: %.2f%%  |  Mann-Whitney p = %.6g\n\n", delta_pct, wt$p.value))

  list(loco = loco_results, random = random_results)
}

## ============================================================
## Run both models
## ============================================================
rf_results  <- run_loco_and_baseline(fit_predict_rf,  "rf",  n_reps_random = 20)
svr_results <- run_loco_and_baseline(fit_predict_svr, "svr", n_reps_random = 20)

## ============================================================
## DIRECT COMPARISON: does spatial adaptation help under genuine
## spatial extrapolation? Paired by city against the SRF/SSVR LOCO
## results you already have on disk.
## ============================================================
srf_path  <- file.path(PROJECT_ROOT, "outputs", "05_srf_spatial_loco", "srf_spatial_loco.csv")
ssvr_path <- file.path(PROJECT_ROOT, "outputs", "06_ssvr_spatial_loco", "ssvr_spatial_loco.csv")

if (file.exists(srf_path) && file.exists(ssvr_path)) {

  srf_loco_raw  <- read.csv(srf_path)
  ssvr_loco_raw <- read.csv(ssvr_path)

  ## FIX: don't assume the four data frames already have the same rows,
  ## in the same order, with no duplicates/NAs -- that assumption is what
  ## caused "'x' and 'y' must have the same length". Instead:
  ##   1. keep only finite-RMSE rows and drop duplicate City entries
  ##   2. explicitly merge all four on City, so only cities present and
  ##      valid in ALL FOUR files are compared (an inner join)
  clean_loco <- function(df) {
    df <- df[is.finite(df$RMSE), ]
    df <- df[!duplicated(df$City), ]
    df
  }

  rf_loco   <- clean_loco(rf_results$loco)
  svr_loco  <- clean_loco(svr_results$loco)
  srf_loco  <- clean_loco(srf_loco_raw)
  ssvr_loco <- clean_loco(ssvr_loco_raw)

  cat(sprintf("\nRows going into the comparison (before aligning): RF=%d SRF=%d SVR=%d SSVR=%d\n",
              nrow(rf_loco), nrow(srf_loco), nrow(svr_loco), nrow(ssvr_loco)))

  ## Inner-join RF+SRF pair and SVR+SSVR pair separately by City, so each
  ## comparison uses exactly the matched, aligned set of cities.
  rf_srf   <- merge(rf_loco[, c("City","RMSE")], srf_loco[, c("City","RMSE")],
                     by = "City", suffixes = c("_RF", "_SRF"))
  svr_ssvr <- merge(svr_loco[, c("City","RMSE")], ssvr_loco[, c("City","RMSE")],
                     by = "City", suffixes = c("_SVR", "_SSVR"))

  cat(sprintf("Cities matched for RF vs SRF: %d | for SVR vs SSVR: %d\n",
              nrow(rf_srf), nrow(svr_ssvr)))
  if (nrow(rf_srf) < 18 || nrow(svr_ssvr) < 18) {
    warning("Fewer than 18 cities matched -- check the *_spatial_loco.csv files for ",
            "missing/duplicate/NA City entries before trusting the test below.")
  }

  cat("\n==================== CORE QUESTION: does spatial adaptation help under LOCO? ====================\n")

  cat(sprintf("\nRF   mean LOCO RMSE = %.4f\nSRF  mean LOCO RMSE = %.4f\n",
              mean(rf_srf$RMSE_RF), mean(rf_srf$RMSE_SRF)))
  w1 <- wilcox.test(rf_srf$RMSE_SRF, rf_srf$RMSE_RF, paired = TRUE, alternative = "less")
  cat(sprintf("Paired Wilcoxon (SRF < RF), one-sided: W=%.2f p=%.5f\n", w1$statistic, w1$p.value))
  cat(sprintf("SRF better than RF in %d of %d cities\n",
              sum(rf_srf$RMSE_SRF < rf_srf$RMSE_RF), nrow(rf_srf)))

  cat(sprintf("\nSVR  mean LOCO RMSE = %.4f\nSSVR mean LOCO RMSE = %.4f\n",
              mean(svr_ssvr$RMSE_SVR), mean(svr_ssvr$RMSE_SSVR)))
  w2 <- wilcox.test(svr_ssvr$RMSE_SSVR, svr_ssvr$RMSE_SVR, paired = TRUE, alternative = "less")
  cat(sprintf("Paired Wilcoxon (SSVR < SVR), one-sided: W=%.2f p=%.5f\n", w2$statistic, w2$p.value))
  cat(sprintf("SSVR better than SVR in %d of %d cities\n",
              sum(svr_ssvr$RMSE_SSVR < svr_ssvr$RMSE_SVR), nrow(svr_ssvr)))

  write.csv(rf_srf,   "rf_vs_srf_loco_comparison.csv",   row.names = FALSE)
  write.csv(svr_ssvr, "svr_vs_ssvr_loco_comparison.csv", row.names = FALSE)
  cat("\nSaved: rf_vs_srf_loco_comparison.csv, svr_vs_ssvr_loco_comparison.csv\n")

} else {
  cat("\nNOTE: srf_spatial_loco.csv and/or ssvr_spatial_loco.csv not found at\n")
  cat(sprintf("  %s\n  %s\n", srf_path, ssvr_path))
  cat("Run 05_srf_spatial_loco and 06_ssvr_spatial_loco first to generate them.\n")
}
