## Set the path to your desired folder with this command
PROJECT_ROOT <- Sys.getenv("IRAQ_WEATHER_ROOT", unset = getwd())  # set by run_all.R; see README
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "08_sensitivity_analysis")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)
getwd()



## ============================================================
## Sensitivity analysis: does the core RF-vs-SRF LOCO finding
## (RF outperforms SRF under leave-one-governorate-out) hold up to
## changes in (1) random seed, (2) hyperparameter search grid, and
## (3) the spatial-weighting scheme used to build the lag predictors?
##
## Design: rather than a full factorial (which would be extremely
## expensive), we vary ONE factor at a time relative to the baseline
## configuration already reported in the paper (Seed=123, narrow grid,
## uniform k-NN weights):
##   - Baseline                  : Seed=123, narrow grid, uniform weights
##   - Seed sensitivity          : Seed=456, Seed=789 (same grid, weights)
##   - Hyperparameter sensitivity: wider ntree/mtry grid (same seed, weights)
##   - Spatial-weight sensitivity: inverse-distance weighting instead of
##                                 uniform k-NN averaging (same seed, grid)
##
## RUNTIME WARNING: each configuration re-runs the full 18-fold LOCO
## comparison with a 12-fold CV grid search per fold. This script runs
## 5 configurations total (1 baseline + 4 variants) -- expect this to
## take several times as long as a single run of the original RF-vs-SRF
## comparison. Checkpointing is included per configuration so an
## interruption does not lose completed configurations.
## ============================================================

library(openxlsx)
library(spdep)
library(randomForest)
library(RANN)

path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
MyData <- read.xlsx(path)
colnames(MyData) <- c("ID","Y","X1","X2","X3","X4","X5","Latitude","Longitude")

Formula        <- Y ~ X1 + X2 + X3 + X4 + X5
SpatialFormula <- Y ~ X1 + X2 + X3 + X4 + X5 + WY + WX1 + WX2 + WX3 + WX4 + WX5

## ============================================================
## build_city_lag(): now parameterized by `weight_scheme`.
##   "uniform" (the original/baseline scheme): each of the k nearest
##     cities gets equal weight 1/k (spdep's default style="W" on an
##     unweighted neighbor list).
##   "idw": each of the k nearest cities is weighted by inverse distance
##     (closer neighboring cities contribute more), then row-standardized.
## ============================================================
build_city_lag <- function(TrainData_base, TestData_base, MyData_full,
                            k_lag_city, weight_scheme = "uniform") {

  train_city_stats <- aggregate(cbind(Y, X1, X2, X3, X4, X5) ~ ID, data = TrainData_base, FUN = mean)
  train_city_coords_lookup <- aggregate(cbind(Longitude, Latitude) ~ ID, data = TrainData_base, FUN = function(x) x[1])
  train_city_stats <- merge(train_city_stats, train_city_coords_lookup, by = "ID")
  train_city_stats <- train_city_stats[order(train_city_stats$ID), ]

  n_cities_train <- nrow(train_city_stats)
  k_eff <- min(k_lag_city, n_cities_train - 1)

  city_coords_mat <- as.matrix(train_city_stats[, c("Longitude","Latitude")])
  city_nb <- knn2nb(knearneigh(city_coords_mat, k = k_eff))

  if (weight_scheme == "uniform") {
    city_weights <- nb2listw(city_nb, style = "W")
  } else if (weight_scheme == "idw") {
    nb_dists <- nbdists(city_nb, city_coords_mat)
    idw_glist <- lapply(nb_dists, function(d) 1 / pmax(d, 1e-6))  # avoid divide-by-zero
    city_weights <- nb2listw(city_nb, glist = idw_glist, style = "W")  # style="W" row-standardizes AFTER applying idw weights
  } else {
    stop("Unknown weight_scheme: ", weight_scheme)
  }

  train_city_stats$WY_city  <- lag.listw(city_weights, train_city_stats$Y)
  train_city_stats$WX1_city <- lag.listw(city_weights, train_city_stats$X1)
  train_city_stats$WX2_city <- lag.listw(city_weights, train_city_stats$X2)
  train_city_stats$WX3_city <- lag.listw(city_weights, train_city_stats$X3)
  train_city_stats$WX4_city <- lag.listw(city_weights, train_city_stats$X4)
  train_city_stats$WX5_city <- lag.listw(city_weights, train_city_stats$X5)

  city_lag_cols <- c("ID","WY_city","WX1_city","WX2_city","WX3_city","WX4_city","WX5_city")
  TrainData_k <- merge(TrainData_base, train_city_stats[, city_lag_cols], by = "ID", all.x = TRUE)
  TestData_k  <- merge(TestData_base,  train_city_stats[, city_lag_cols], by = "ID", all.x = TRUE)

  missing_test_cities <- unique(TestData_k$ID[is.na(TestData_k$WY_city)])
  if (length(missing_test_cities) > 0) {
    all_city_coords <- aggregate(cbind(Longitude, Latitude) ~ ID, data = MyData_full, FUN = function(x) x[1])
    for (missing_id in missing_test_cities) {
      this_coord <- as.matrix(all_city_coords[all_city_coords$ID == missing_id, c("Longitude","Latitude")])
      nn <- RANN::nn2(data = city_coords_mat, query = this_coord, k = 1)
      nearest_row <- train_city_stats[nn$nn.idx[1], ]
      fill_idx <- TestData_k$ID == missing_id
      TestData_k$WY_city[fill_idx]  <- nearest_row$WY_city
      TestData_k$WX1_city[fill_idx] <- nearest_row$WX1_city
      TestData_k$WX2_city[fill_idx] <- nearest_row$WX2_city
      TestData_k$WX3_city[fill_idx] <- nearest_row$WX3_city
      TestData_k$WX4_city[fill_idx] <- nearest_row$WX4_city
      TestData_k$WX5_city[fill_idx] <- nearest_row$WX5_city
    }
  }

  TrainData_k$WY  <- TrainData_k$WY_city;  TrainData_k$WX1 <- TrainData_k$WX1_city
  TrainData_k$WX2 <- TrainData_k$WX2_city; TrainData_k$WX3 <- TrainData_k$WX3_city
  TrainData_k$WX4 <- TrainData_k$WX4_city; TrainData_k$WX5 <- TrainData_k$WX5_city
  TestData_k$WY  <- TestData_k$WY_city;  TestData_k$WX1 <- TestData_k$WX1_city
  TestData_k$WX2 <- TestData_k$WX2_city; TestData_k$WX3 <- TestData_k$WX3_city
  TestData_k$WX4 <- TestData_k$WX4_city; TestData_k$WX5 <- TestData_k$WX5_city

  list(TrainData = TrainData_k, TestData = TestData_k)
}

## ============================================================
## fit_predict_rf() and fit_predict_srf(): parameterized by seed and
## hyperparameter grid, so the same functions serve both the baseline
## and every sensitivity variant.
## ============================================================
fit_predict_rf <- function(TrainData, TestData, seed, ntree.grid, mtry.grid) {
  set.seed(seed)
  K <- 12
  n_k <- nrow(TrainData)
  fold.id <- sample(rep(1:K, length.out = n_k))
  grid.results <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())
  for (nt in ntree.grid) for (mt in mtry.grid) {
    rmse.fold <- numeric(K)
    for (kk in 1:K) {
      Train.CV <- TrainData[fold.id != kk, ]; Valid.CV <- TrainData[fold.id == kk, ]
      set.seed(seed)
      rf.cv <- randomForest(Formula, data = Train.CV, ntree = nt, mtry = mt)
      pred.cv <- predict(rf.cv, newdata = Valid.CV)
      rmse.fold[kk] <- sqrt(mean((Valid.CV$Y - pred.cv)^2))
    }
    grid.results <- rbind(grid.results, data.frame(ntree = nt, mtry = mt, RMSE = mean(rmse.fold)))
  }
  best <- grid.results[which.min(grid.results$RMSE), ]
  model <- randomForest(Formula, data = TrainData, ntree = best$ntree, mtry = best$mtry)
  pred <- predict(model, newdata = TestData)
  sqrt(mean((TestData$Y - pred)^2))
}

fit_predict_srf <- function(TrainData, TestData, seed, ntree.grid, mtry.grid,
                             k_lag_city, weight_scheme, MyData_full) {
  lag_result <- build_city_lag(TrainData, TestData, MyData_full, k_lag_city, weight_scheme)
  TrainData_k <- lag_result$TrainData; TestData_k <- lag_result$TestData
  set.seed(seed)
  K <- 12
  n_k <- nrow(TrainData_k)
  fold.id <- sample(rep(1:K, length.out = n_k))
  grid.results <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())
  for (nt in ntree.grid) for (mt in mtry.grid) {
    rmse.fold <- numeric(K)
    for (kk in 1:K) {
      Train.CV <- TrainData_k[fold.id != kk, ]; Valid.CV <- TrainData_k[fold.id == kk, ]
      set.seed(seed)
      srf.cv <- randomForest(SpatialFormula, data = Train.CV, ntree = nt, mtry = mt)
      pred.cv <- predict(srf.cv, newdata = Valid.CV)
      rmse.fold[kk] <- sqrt(mean((Valid.CV$Y - pred.cv)^2))
    }
    grid.results <- rbind(grid.results, data.frame(ntree = nt, mtry = mt, RMSE = mean(rmse.fold)))
  }
  best <- grid.results[which.min(grid.results$RMSE), ]
  model <- randomForest(SpatialFormula, data = TrainData_k, ntree = best$ntree, mtry = best$mtry)
  pred <- predict(model, newdata = TestData_k)
  sqrt(mean((TestData_k$Y - pred)^2))
}

## ============================================================
## run_loco_comparison(): the full 18-fold RF-vs-SRF LOCO comparison
## for one configuration, with checkpointing keyed by `config_label`.
## ============================================================
run_loco_comparison <- function(config_label, seed, ntree.grid, mtry.grid,
                                 k_lag_city = 5, weight_scheme = "uniform") {

  checkpoint_file <- paste0("sensitivity_", config_label, ".csv")
  all_cities <- sort(unique(MyData$ID))

  already_done <- integer(0)
  if (file.exists(checkpoint_file)) {
    prev <- read.csv(checkpoint_file)
    already_done <- prev$City[is.finite(prev$RF) & is.finite(prev$SRF)]
    cat(sprintf("[%s] Resuming: %d/%d cities already done\n", config_label, length(already_done), length(all_cities)))
  } else {
    write.csv(data.frame(City = integer(0), RF = numeric(0), SRF = numeric(0)),
              checkpoint_file, row.names = FALSE)
  }

  for (city_id in all_cities) {
    if (city_id %in% already_done) next
    TrainData_fold <- MyData[MyData$ID != city_id, ]
    TestData_fold  <- MyData[MyData$ID == city_id, ]

    rf_rmse  <- tryCatch(fit_predict_rf(TrainData_fold, TestData_fold, seed, ntree.grid, mtry.grid),
                          error = function(e) NA_real_)
    srf_rmse <- tryCatch(fit_predict_srf(TrainData_fold, TestData_fold, seed, ntree.grid, mtry.grid,
                                         k_lag_city, weight_scheme, MyData),
                          error = function(e) NA_real_)

    cat(sprintf("[%s] city=%d: RF=%.4f SRF=%.4f\n", config_label, city_id, rf_rmse, srf_rmse))
    write.table(data.frame(City = city_id, RF = rf_rmse, SRF = srf_rmse),
                checkpoint_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }

  results <- read.csv(checkpoint_file)
  results <- results[!duplicated(results$City), ]
  results <- results[is.finite(results$RF) & is.finite(results$SRF), ]
  write.csv(results, checkpoint_file, row.names = FALSE)

  wt <- wilcox.test(results$SRF, results$RF, alternative = "less")  # does SRF < RF (SRF wins)?
  diffs <- results$SRF - results$RF
  ranks <- rank(abs(diffs))
  r_rb <- (sum(ranks[diffs > 0]) - sum(ranks[diffs < 0])) / sum(ranks)

  cat(sprintf("\n=== %s: RF mean=%.4f, SRF mean=%.4f, W=%.1f, p=%.4f, r=%.3f, SRF-better-in=%d/%d ===\n\n",
              config_label, mean(results$RF), mean(results$SRF), wt$statistic, wt$p.value, r_rb,
              sum(diffs < 0), nrow(results)))

  data.frame(config = config_label, RF_mean = mean(results$RF), SRF_mean = mean(results$SRF),
             W = wt$statistic, p = wt$p.value, r = r_rb, n_srf_better = sum(diffs < 0), n = nrow(results))
}

## ============================================================
## Run baseline + 4 sensitivity variants
## ============================================================
NARROW_NTREE <- c(100, 500, 1000); NARROW_MTRY <- c(2, 3, 4)
WIDE_NTREE   <- c(100, 300, 500, 800, 1000, 1500); WIDE_MTRY <- c(2, 3, 4, 5, 6)

summary_all <- data.frame()

summary_all <- rbind(summary_all, run_loco_comparison(
  "baseline_seed123", seed = 123, ntree.grid = NARROW_NTREE, mtry.grid = NARROW_MTRY))

summary_all <- rbind(summary_all, run_loco_comparison(
  "seed456", seed = 456, ntree.grid = NARROW_NTREE, mtry.grid = NARROW_MTRY))

summary_all <- rbind(summary_all, run_loco_comparison(
  "seed789", seed = 789, ntree.grid = NARROW_NTREE, mtry.grid = NARROW_MTRY))

summary_all <- rbind(summary_all, run_loco_comparison(
  "wide_grid", seed = 123, ntree.grid = WIDE_NTREE, mtry.grid = WIDE_MTRY))

summary_all <- rbind(summary_all, run_loco_comparison(
  "idw_weights", seed = 123, ntree.grid = NARROW_NTREE, mtry.grid = NARROW_MTRY,
  weight_scheme = "idw"))

cat("\n\n==================== FULL SENSITIVITY SUMMARY ====================\n")
print(summary_all, row.names = FALSE)
write.csv(summary_all, "sensitivity_full_summary.csv", row.names = FALSE)
cat("\nSaved: sensitivity_full_summary.csv\n")

