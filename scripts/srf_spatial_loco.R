## ============================================================
## Spatial LOCO validation for SRF, mirroring the SSVR analysis just
## completed, using the leakage-free city-level spatial lag pipeline.
##
## RUNTIME WARNING (READ BEFORE RUNNING): SRF's hyperparameter search is
## a full grid search (3 ntree x 3 mtry = 9 combos) with a MANUAL 12-fold
## CV loop per combo -- i.e. 108 randomForest fits PER FOLD, vs. SSVR's
## much cheaper svm tuning. With 18 LOCO folds + N_REPS_RANDOM random
## splits, that is roughly (18 + N_REPS_RANDOM) x 108 randomForest fits
## in total. This will likely take substantially longer than the SSVR
## version. Suggestions if this is too slow:
##   - reduce N_REPS_RANDOM below (default lowered to 10, not 20, for
##     this reason)
##   - reduce ntree.grid to c(100, 500) (drops 1000-tree fits)
##   - run overnight / unattended; checkpointing means an interruption
##     does not lose completed folds
## ============================================================

PROJECT_ROOT <- Sys.getenv("IRAQ_WEATHER_ROOT", unset = getwd())  # set by run_all.R; see README
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "05_srf_spatial_loco")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)

library(openxlsx)
library(sf)
library(spdep)
library(randomForest)
library(RANN)

Seed <- 123
K_LAG_CITY <- 5   # same default used for the main SRF result in Table 2

path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
MyData <- read.xlsx(path)
colnames(MyData) <- c("ID","Y","X1","X2","X3","X4","X5","Latitude","Longitude")

## ============================================================
## build_city_lag(): identical to the function in the main pipeline and
## in ssvr_spatial_loco.R. Included here so this script is self-contained.
## ============================================================
build_city_lag <- function(TrainData_base, TestData_base, MyData_full, k_lag_city) {

  train_city_stats = aggregate(cbind(Y, X1, X2, X3, X4, X5) ~ ID, data = TrainData_base, FUN = mean)
  train_city_coords_lookup = aggregate(cbind(Longitude, Latitude) ~ ID, data = TrainData_base, FUN = function(x) x[1])
  train_city_stats = merge(train_city_stats, train_city_coords_lookup, by = "ID")
  train_city_stats = train_city_stats[order(train_city_stats$ID), ]

  n_cities_train = nrow(train_city_stats)
  k_eff = min(k_lag_city, n_cities_train - 1)

  city_coords_mat = as.matrix(train_city_stats[, c("Longitude","Latitude")])
  city_nb      = knn2nb(knearneigh(city_coords_mat, k = k_eff))
  city_weights = nb2listw(city_nb, style = "W")

  train_city_stats$WY_city  = lag.listw(city_weights, train_city_stats$Y)
  train_city_stats$WX1_city = lag.listw(city_weights, train_city_stats$X1)
  train_city_stats$WX2_city = lag.listw(city_weights, train_city_stats$X2)
  train_city_stats$WX3_city = lag.listw(city_weights, train_city_stats$X3)
  train_city_stats$WX4_city = lag.listw(city_weights, train_city_stats$X4)
  train_city_stats$WX5_city = lag.listw(city_weights, train_city_stats$X5)

  city_lag_cols = c("ID","WY_city","WX1_city","WX2_city","WX3_city","WX4_city","WX5_city")
  TrainData_k = merge(TrainData_base, train_city_stats[, city_lag_cols], by = "ID", all.x = TRUE)
  TestData_k  = merge(TestData_base,  train_city_stats[, city_lag_cols], by = "ID", all.x = TRUE)

  missing_test_cities = unique(TestData_k$ID[is.na(TestData_k$WY_city)])
  if (length(missing_test_cities) > 0) {
    all_city_coords = aggregate(cbind(Longitude, Latitude) ~ ID, data = MyData_full, FUN = function(x) x[1])
    for (missing_id in missing_test_cities) {
      this_coord = as.matrix(all_city_coords[all_city_coords$ID == missing_id, c("Longitude","Latitude")])
      nn = RANN::nn2(data = city_coords_mat, query = this_coord, k = 1)
      nearest_row = train_city_stats[nn$nn.idx[1], ]
      fill_idx = TestData_k$ID == missing_id
      TestData_k$WY_city[fill_idx]  = nearest_row$WY_city
      TestData_k$WX1_city[fill_idx] = nearest_row$WX1_city
      TestData_k$WX2_city[fill_idx] = nearest_row$WX2_city
      TestData_k$WX3_city[fill_idx] = nearest_row$WX3_city
      TestData_k$WX4_city[fill_idx] = nearest_row$WX4_city
      TestData_k$WX5_city[fill_idx] = nearest_row$WX5_city
    }
  }

  TrainData_k$WY  = TrainData_k$WY_city;  TrainData_k$WX1 = TrainData_k$WX1_city
  TrainData_k$WX2 = TrainData_k$WX2_city; TrainData_k$WX3 = TrainData_k$WX3_city
  TrainData_k$WX4 = TrainData_k$WX4_city; TrainData_k$WX5 = TrainData_k$WX5_city
  TestData_k$WY  = TestData_k$WY_city;  TestData_k$WX1 = TestData_k$WX1_city
  TestData_k$WX2 = TestData_k$WX2_city; TestData_k$WX3 = TestData_k$WX3_city
  TestData_k$WX4 = TestData_k$WX4_city; TestData_k$WX5 = TestData_k$WX5_city

  list(TrainData = TrainData_k, TestData = TestData_k)
}

SpatialFormula = Y ~ X1 + X2 + X3 + X4 + X5 + WY + WX1 + WX2 + WX3 + WX4 + WX5

## ---- Full 12-fold CV grid search + fit + predict, for one train/test
## partition (used identically for both LOCO folds and random-split reps) ----
fit_predict_srf <- function(TrainData_base, TestData_base, label,
                             ntree.grid = c(100, 500, 1000), mtry.grid = c(2, 3, 4)) {

  lag_result <- build_city_lag(TrainData_base, TestData_base, MyData, K_LAG_CITY)
  TrainData_k <- lag_result$TrainData
  TestData_k  <- lag_result$TestData

  set.seed(Seed)
  K <- 12
  n_k <- nrow(TrainData_k)
  fold.id <- sample(rep(1:K, length.out = n_k))

  grid.results <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())
  for (nt in ntree.grid) {
    for (mt in mtry.grid) {
      rmse.fold <- numeric(K)
      for (kk in 1:K) {
        Train.CV <- TrainData_k[fold.id != kk, ]
        Valid.CV <- TrainData_k[fold.id == kk, ]
        set.seed(Seed)
        rf.cv <- randomForest(SpatialFormula, data = Train.CV, ntree = nt, mtry = mt)
        pred.cv <- predict(rf.cv, newdata = Valid.CV)
        rmse.fold[kk] <- sqrt(mean((Valid.CV$Y - pred.cv)^2))
      }
      grid.results <- rbind(grid.results, data.frame(ntree = nt, mtry = mt, RMSE = mean(rmse.fold)))
    }
  }
  best <- grid.results[which.min(grid.results$RMSE), ]

  model <- randomForest(SpatialFormula, data = TrainData_k, ntree = best$ntree, mtry = best$mtry, importance = TRUE)
  pred  <- predict(model, newdata = TestData_k)

  rmse <- sqrt(mean((TestData_k$Y - pred)^2))
  mae  <- mean(abs(TestData_k$Y - pred))
  r2   <- 1 - sum((TestData_k$Y - pred)^2) / sum((TestData_k$Y - mean(TestData_k$Y))^2)

  cat(sprintf("[%s] R2=%.4f RMSE=%.4f MAE=%.4f (ntree=%d, mtry=%d, n_train=%d, n_test=%d)\n",
              label, r2, rmse, mae, best$ntree, best$mtry, nrow(TrainData_k), nrow(TestData_k)))

  list(R2 = r2, RMSE = rmse, MAE = mae)
}

## ============================================================
## PART 1: Spatial LOCO -- hold out each of the 18 cities in turn.
## Checkpointed: safe to interrupt and resume.
## ============================================================
all_cities <- sort(unique(MyData$ID))
checkpoint_loco <- "srf_spatial_loco.csv"

already_done <- integer(0)
if (file.exists(checkpoint_loco)) {
  prev <- read.csv(checkpoint_loco)
  already_done <- prev$City[is.finite(prev$RMSE)]
  cat(sprintf("Resuming: %d cities already done, %d remaining\n",
              length(already_done), length(all_cities) - length(already_done)))
} else {
  write.csv(data.frame(City=integer(0), R2=numeric(0), RMSE=numeric(0), MAE=numeric(0)),
            checkpoint_loco, row.names = FALSE)
}

for (city_id in all_cities) {
  if (city_id %in% already_done) next

  TrainData_fold <- MyData[MyData$ID != city_id, ]
  TestData_fold  <- MyData[MyData$ID == city_id, ]

  res <- tryCatch(
    fit_predict_srf(TrainData_fold, TestData_fold, sprintf("SRF-LOCO city=%d", city_id)),
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

cat("\n==================== SRF SPATIAL LOCO SUMMARY ====================\n")
print(loco_results, row.names = FALSE)
cat(sprintf("\nMean RMSE = %.4f, SD = %.4f (n=%d cities)\n",
            mean(loco_results$RMSE, na.rm=TRUE), sd(loco_results$RMSE, na.rm=TRUE),
            sum(is.finite(loco_results$RMSE))))

## ============================================================
## PART 2: Repeated random-split baseline for SRF (same k=5 city-level
## lag construction, same 75/25 ratio as Table 2).
##
## N_REPS_RANDOM lowered to 10 (not 20, as used for SSVR) given SRF's
## higher per-fold cost -- raise it back to 20 if you have the time and
## want a more precise sd estimate.
## ============================================================
N_REPS_RANDOM <- 10
checkpoint_random <- "srf_random_split_baseline.csv"

already_done_r <- integer(0)
if (file.exists(checkpoint_random)) {
  prev_r <- read.csv(checkpoint_random)
  already_done_r <- prev_r$rep_id[is.finite(prev_r$RMSE)]
  cat(sprintf("Resuming: %d random-split reps already done\n", length(already_done_r)))
} else {
  write.csv(data.frame(rep_id=integer(0), R2=numeric(0), RMSE=numeric(0), MAE=numeric(0)),
            checkpoint_random, row.names = FALSE)
}

n_total <- nrow(MyData)
for (r in seq_len(N_REPS_RANDOM)) {
  if (r %in% already_done_r) next

  set.seed(20240715 + r)
  idx <- sample(1:n_total, size = round(0.75 * n_total))
  TrainData_r <- MyData[idx, ]
  TestData_r  <- MyData[-idx, ]

  res_r <- tryCatch(
    fit_predict_srf(TrainData_r, TestData_r, sprintf("SRF-Random rep=%d", r)),
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

cat("\n==================== SRF RANDOM-SPLIT BASELINE SUMMARY ====================\n")
cat(sprintf("Mean RMSE = %.4f, SD = %.4f (n=%d reps)\n",
            mean(random_results$RMSE, na.rm=TRUE), sd(random_results$RMSE, na.rm=TRUE),
            sum(is.finite(random_results$RMSE))))

## ============================================================
## PART 3: Formal comparison, mirroring the SSVR analysis
## ============================================================
wt <- wilcox.test(loco_results$RMSE, random_results$RMSE, alternative = "greater")
cat(sprintf("\nMann-Whitney (SRF Spatial LOCO > SRF Random split): W=%s, p=%.6g\n",
            format(wt$statistic), wt$p.value))

delta_pct <- (mean(loco_results$RMSE, na.rm=TRUE) - mean(random_results$RMSE, na.rm=TRUE)) /
             mean(random_results$RMSE, na.rm=TRUE) * 100
cat(sprintf("Relative RMSE increase (LOCO vs random): %.2f%%\n", delta_pct))
