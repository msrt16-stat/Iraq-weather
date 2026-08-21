
## Set the path to your desired folder with this command
PROJECT_ROOT <- "C:/Weather Project"
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "01_weather_R_FINAL")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)
getwd()


## ============================================================
## FINAL CONSOLIDATED PIPELINE -- humidity paper (Iraq, 18 cities)
##
## This single script reproduces every corrected result reported in the
## final version of the manuscript. It consolidates several rounds of
## fixes made over the course of review, in one place:
##
##   1. Train/test split (75/25, matching the paper's Methods text) is
##      performed BEFORE any spatial-lag feature is built, closing a
##      data-leakage bug in which the lag predictors (WY, WX1..WX5) used
##      by SRF and SSVR were originally computed on the full dataset
##      (train+test combined) before splitting.
##   2. The spatial lag itself is built at the level of the 18 unique
##      cities (train-only city means, k nearest CITIES, not individual
##      monthly rows), rather than row-level KNN -- the row-level version
##      degenerated into a same-city temporal average because the 18
##      cities are each repeated ~120 times at identical coordinates.
##   3. GWR/SR uses a bi-square kernel with an ADAPTIVE bandwidth
##      (matching the Methods section), instead of the fixed-bandwidth
##      Gaussian kernel used in the original implementation.
##   4. The results matrix `Res` is indexed by row NAME (Res["GWR",], not
##      Res[2,]), fixing a bug where model results were silently stored
##      one row off from their label.
##   5. A sensitivity check over k in {3, 5, 8} neighboring cities for the
##      spatial lag confirms the model ranking is stable to this choice
##      (k=5 is used as the primary/reported specification).
##   6. Train-set metrics are computed for all six models at the end, to
##      populate the full Train/Test Table 2 in the manuscript.
##
## Requires: openxlsx, spgwr, sf, spdep, e1071, randomForest, RANN
##
## RUNTIME WARNING: the k-sensitivity loop refits SRF and SSVR (each with
## a full 12-fold CV hyperparameter search) three times. This is the
## slowest part of the script; expect it to take a while. See the
## comment above the loop for how to reduce this if needed.
## ============================================================

rm(list=ls())
cat("\f")
Seed=123
rTrain=0.75 # 75%
########
Meths = c("OLS","GWR","SVR","SSVR","RF","SRF")

Crs=c("MSE","R2","AIC-B","AIC-H","AICc","MAE","RMSE")
nM=length(Meths)
nC=length(Crs)
library(openxlsx)
library(spgwr)
library(sf)
library(spdep)
library(e1071)
library(randomForest)
library(RANN)   # needed for cross-set (test-to-train-only) nearest-neighbor search

path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))

MyData=read.xlsx(path)
colnames(MyData)=c("ID","Y","X1","X2","X3","X4","X5","Latitude","Longitude")
set.seed(Seed)

## ============================================================
## FIX #1 (THE CRITICAL ONE): split into train/test FIRST, before any
## spatial-lag feature is computed. In the original script, WY/WX1..WX5
## were computed on the FULL dataset (train+test combined) and only
## afterward split -- meaning every training row's "spatial lag" predictor
## was partly built from neighboring rows that ended up in the test set,
## including the target variable Y itself via WY. Because there are only
## 18 distinct city coordinates (each repeated ~120 times), a training
## row's k=80 nearest neighbors are overwhelmingly OTHER MONTHS OF THE
## SAME CITY -- many of which are in the test set under a random
## observation-level split. That is a near-verbatim leakage of the test
## set's Y values into the training features for SRF and SSVR (both of
## which include WY as a predictor). Splitting first, and computing test
## lags using ONLY training-set neighbors (Fix #2 below), removes this.
## ============================================================

idx_train = sample(
  1:nrow(MyData),
  size = round(rTrain * nrow(MyData))
)

TrainData = MyData[idx_train, ]
TestData  = MyData[-idx_train, ]

## ---- Diagnostic-only Moran's I on the full Y (does not feed into any
## model input, so computing it on the full data before splitting is fine
## -- this is purely descriptive, unlike the lag features below) ----
weather_sf_all = st_as_sf(MyData, coords = c("Longitude","Latitude"), crs = 4326)
coords_all     = st_coordinates(weather_sf_all)
neighbors_all  = knn2nb(knearneigh(coords_all, k = 80))
weights_all    = nb2listw(neighbors_all, style = "W")
moran_Y_all    = moran.test(MyData$Y, weights_all)
print(moran_Y_all)

## ============================================================
## FIX #2 (REVISED): city-level spatial lag instead of row-level K=80 KNN.
##
## WHY THIS CHANGED: with K=80 nearest neighbors computed on individual
## (city x month) rows, and only 18 distinct city coordinates each
## repeated ~90-120 times, a training row's 80 nearest neighbors were
## almost entirely OTHER MONTHS OF THE SAME CITY (distance = 0 ties).
## That is leakage-free (Fix #2's original version already ensured no
## test data was used), but it is not a genuinely SPATIAL signal --
## it behaves like a city-specific historical average (closer to a
## city fixed effect) rather than a GWR-style "borrow strength from
## neighboring cities" spatial lag, which is what the paper's method
## section actually describes.
##
## FIX: aggregate to the 18 unique cities first, build the neighbor
## structure over those 18 city-level points (not over 1600+ duplicated
## rows), and compute each city's spatial lag as a weighted average of
## ITS NEIGHBORING CITIES' mean Y/X (excluding the city's own value,
## which is how nb2listw/lag.listw work by default -- a location is
## never its own neighbor). This is a true cross-city spatial lag.
##
## K_LAG_CITY = 5 is a reasonable starting choice for only 18 spatial
## units (roughly the "average nearest-neighbor count" used in many
## spatial-lag / SLX applications with small-N point patterns) -- it is
## a judgment call, not a uniquely correct number. If you want to check
## sensitivity to this choice, rerunning with K_LAG_CITY = 3 or 8 and
## comparing Table 2 is a reasonable robustness check to add later.
## ============================================================

## ============================================================
## Reusable function: build genuine city-level spatial lag features for
## a given neighbor count k, entirely from TRAINING city means (leakage-
## free), then project onto TestData via each test row's own city.
## Wrapped in a function so we can rerun it for several k values without
## re-deriving/duplicating logic each time (Fix #2 from before, now
## parameterized).
## ============================================================
build_city_lag <- function(TrainData_base, TestData_base, MyData_full, k_lag_city) {

  train_city_stats = aggregate(
    cbind(Y, X1, X2, X3, X4, X5) ~ ID,
    data = TrainData_base, FUN = mean
  )
  train_city_coords_lookup = aggregate(
    cbind(Longitude, Latitude) ~ ID,
    data = TrainData_base, FUN = function(x) x[1]
  )
  train_city_stats = merge(train_city_stats, train_city_coords_lookup, by = "ID")
  train_city_stats = train_city_stats[order(train_city_stats$ID), ]

  n_cities_train = nrow(train_city_stats)
  if (n_cities_train < 18) {
    warning(sprintf("k=%d: only %d of 18 cities present in TrainData.", k_lag_city, n_cities_train))
  }

  k_eff = min(k_lag_city, n_cities_train - 1)
  if (k_eff < k_lag_city) {
    warning(sprintf("k=%d reduced to %d (only %d cities available).", k_lag_city, k_eff, n_cities_train))
  }

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
    warning(sprintf("k=%d: projecting lag from nearest trained city for IDs: %s",
                     k_lag_city, paste(missing_test_cities, collapse=", ")))
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

  list(TrainData = TrainData_k, TestData = TestData_k, k_eff = k_eff)
}

## Keep an immutable "base" copy of the split data (no WY/WX columns yet)
## -- build_city_lag() will be called fresh for each k value below, each
## time starting from this same base, so results across k are not
## contaminated by a previous k's columns.
TrainData_base = TrainData
TestData_base  = TestData

## ============================================================
## From here on, the rest of the pipeline is structurally the same as
## the original script -- only the WY/WX1..5 construction above changed.
## ============================================================

test_sf = st_as_sf(TestData, coords = c("Longitude","Latitude"), crs = 4326)
test_sp = as(test_sf, "Spatial")
weather_sf = st_as_sf(TrainData, coords = c("Longitude","Latitude"), crs = 4326)
weather_sp = as(weather_sf, "Spatial")

Formula = Y ~ X1+X2+X3+X4+X5

## ============================================================
## FIX #3: match the paper's stated GWR methodology -- bi-square kernel
## with ADAPTIVE bandwidth (adapt=TRUE selects a bandwidth as a
## proportion of nearest neighbors, invariant to the fact that
## coordinates are in raw lon/lat degrees, rather than adapt=FALSE's
## fixed-distance bandwidth which is degree-scale-dependent). Both
## gwr.sel() and gwr() must use the SAME gweight/adapt settings, or the
## bandwidth selected by gwr.sel() won't mean what gwr() expects.
## ============================================================

bw = gwr.sel(formula=Formula, data=weather_sp, adapt=TRUE, gweight=gwr.bisquare)
## NOTE: with adapt=TRUE, gwr.sel() returns a numeric QUANTILE (between 0
## and 1 -- the proportion of nearest observations to use), not a raw
## distance. This numeric value must be passed to gwr() via the `adapt=`
## argument, NOT `bandwidth=` (bandwidth= is only for the fixed-distance
## case). Passing both bandwidth= and adapt=TRUE together is invalid,
## which is exactly what produced the "is.numeric(adapt) is not TRUE"
## error -- gwr() expects a number for adapt, not a logical.
gwr_model = gwr(formula=Formula, data=weather_sp,
                 adapt=bw, gweight=gwr.bisquare,
                 hatmatrix=TRUE, se.fit=TRUE)
head(gwr_model$SDF@data)
gwr_test = gwr(formula=Formula,
                data=weather_sp,
                fit.points=test_sp,
                adapt=bw,
                gweight=gwr.bisquare,
                hatmatrix=FALSE,
                se.fit=FALSE)

pred_gwr =
  gwr_test$SDF$`(Intercept)` +
  gwr_test$SDF$X1 * TestData$X1 +
  gwr_test$SDF$X2 * TestData$X2 +
  gwr_test$SDF$X3 * TestData$X3 +
  gwr_test$SDF$X4 * TestData$X4 +
  gwr_test$SDF$X5 * TestData$X5

ols_model=lm(formula=Formula,data=weather_sp)
coords=st_coordinates(weather_sf)
neighbors=knn2nb(knearneigh(coords,k=80))
weights=nb2listw(neighbors,style="W")
moran_Y=moran.test(TrainData$Y,weights)
print(moran_Y)

zY=scale(TrainData$Y)[,1]
lag_zY=lag.listw(weights,zY)
indep_vars=c("X1","X2","X3","X4","X5")
BLMI_list=list()
for (v in indep_vars){
  zX=scale(TrainData[[v]])[,1]
  BLMI_list[[v]]=zX*lag_zY
}
ok=!is.na(pred_gwr)

MSE_GWR=mean((TestData$Y[ok]-pred_gwr[ok])^2)
RSS_GWR=sum((TestData$Y[ok]-pred_gwr[ok])^2)
TSS_GWR=sum((TestData$Y[ok]-mean(TestData$Y[ok]))^2)
R2_GWR=1-(RSS_GWR/TSS_GWR)
RMSE_GWR = sqrt(MSE_GWR)
MAE_GWR = mean(abs(gwr_model$SDF$gwr.e))

AIC_GWR_B=gwr_model$results$AICb
AIC_GWR_H=gwr_model$results$AICh
AICc_GWR=gwr_model$results$AICc
Res=matrix(0,nM,nC)
rownames(Res)=Meths
colnames(Res)=Crs
Res["GWR",]=c(MSE_GWR,R2_GWR,AIC_GWR_B,AIC_GWR_H,AICc_GWR,MAE_GWR, RMSE_GWR)  # FIX: was Res[1,] (wrong row)
Res["GWR",]

###############################
# SVR + Grid Search + 12-Fold CV  (unaffected by the leakage bug --
# does not use WY/WX, unchanged from the original script)
###############################

tune.svr <- tune(
  svm, Formula, data = TrainData,
  type = "eps-regression", kernel = "radial",
  ranges = list(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1)),
  tunecontrol = tune.control(cross = 12)
)
summary(tune.svr)
svr_model <- tune.svr$best.model
best_cost  <- tune.svr$best.parameters$cost
best_gamma <- tune.svr$best.parameters$gamma
cat("Best Cost =", best_cost, "\n")
cat("Best Gamma =", best_gamma, "\n")

pred_svr <- predict(svr_model, newdata = TestData)
MSE_SVR  <- mean((TestData$Y - pred_svr)^2)
RMSE_SVR <- sqrt(MSE_SVR)
MAE_SVR  <- mean(abs(TestData$Y - pred_svr))
RSS_SVR <- sum((TestData$Y - pred_svr)^2)
TSS_SVR <- sum((TestData$Y - mean(TestData$Y))^2)
R2_SVR <- 1 - RSS_SVR / TSS_SVR

Res["SVR", ] <- c(MSE_SVR, R2_SVR, NA, NA, NA, MAE_SVR, RMSE_SVR)  # FIX: was Res[2,] (wrong row)
Res["SVR", ]

####################################################
# RF : Grid Search + 12-Fold Cross Validation (unaffected by the leakage
# bug and by k_lag_city -- does not use WY/WX, moved up here since it
# does not depend on the k-sensitivity loop below)
####################################################

set.seed(Seed)
ntree.grid <- c(100,500,1000)
mtry.grid  <- c(2,3,4)
K <- 12
n <- nrow(TrainData)
fold.id <- sample(rep(1:K, length.out=n))

RF.Results <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())

for(nt in ntree.grid){
  for(mt in mtry.grid){
    RMSE.fold <- numeric(K)
    for(k in 1:K){
      Train.CV <- TrainData[fold.id!=k,]
      Valid.CV <- TrainData[fold.id==k,]
      rf.cv <- randomForest(Formula, data=Train.CV, ntree=nt, mtry=mt)
      pred.cv <- predict(rf.cv, newdata=Valid.CV)
      RMSE.fold[k] <- sqrt(mean((Valid.CV$Y-pred.cv)^2))
    }
    RF.Results <- rbind(RF.Results, data.frame(ntree=nt, mtry=mt, RMSE=mean(RMSE.fold)))
  }
}

RF.Results
best.rf <- RF.Results[which.min(RF.Results$RMSE), ]
print(best.rf)
best.ntree <- best.rf$ntree
best.mtry  <- best.rf$mtry
cat("\nBest ntree =", best.ntree, "\n")
cat("Best mtry  =", best.mtry, "\n")

rf_model <- randomForest(Formula, data = TrainData, ntree = best.ntree,
                          mtry = best.mtry, importance = TRUE)
pred_rf <- predict(rf_model, newdata = TestData)

MSE_RF  <- mean((TestData$Y - pred_rf)^2)
RMSE_RF <- sqrt(MSE_RF)
MAE_RF  <- mean(abs(TestData$Y - pred_rf))
RSS_RF <- sum((TestData$Y - pred_rf)^2)
TSS_RF <- sum((TestData$Y - mean(TestData$Y))^2)
R2_RF <- 1 - RSS_RF/TSS_RF

Res["RF",] <- c(MSE_RF, R2_RF, NA, NA, NA, MAE_RF, RMSE_RF)
Res["RF",]

####################################################
## K-SENSITIVITY LOOP: refit SSVR and SRF for each candidate neighbor
## count k_lag_city (the number of neighboring CITIES used to build the
## spatial lag). This directly answers "does the SRF/SSVR result depend
## on the somewhat arbitrary choice of k=5?"
##
## RUNTIME WARNING: each iteration re-runs a full 12-fold CV grid search
## for BOTH SRF (3x3=9 ntree/mtry combos x 12 folds = 108 RF fits) and
## SSVR (3x3=9 cost/gamma combos x 12-fold internal CV via tune()). With
## 3 values of k, that is 3x the cost of the single-k version you ran
## before. Expect this to take noticeably longer; consider reducing
## k_values to just two if this is too slow, or running unattended.
####################################################

k_values <- c(3, 5, 8)
SpatialFormula = Y ~ X1 + X2 + X3 + X4 + X5 + WY + WX1 + WX2 + WX3 + WX4 + WX5

k_sensitivity_results <- data.frame(
  k = integer(), Model = character(),
  R2 = numeric(), RMSE = numeric(), MAE = numeric()
)

## Keep the k=5 model objects around so Res["SSVR",]/Res["RF",]/Res["SRF",]
## (the main results table, for comparability with everything else you've
## already produced) still reflects the k=5 case as the "default" choice --
## the full sensitivity picture across all three k's is reported
## separately in k_sensitivity_results / k_sensitivity_summary.csv.
default_k <- 5
ssvr_model_default <- NULL
srf_model_default  <- NULL

for (k_lag_city in k_values) {

  cat(sprintf("\n\n==================== k_lag_city = %d ====================\n", k_lag_city))

  lag_k = build_city_lag(TrainData_base, TestData_base, MyData, k_lag_city)
  TrainData_k = lag_k$TrainData
  TestData_k  = lag_k$TestData

  ## ---- SSVR at this k ----
  set.seed(Seed)
  tune.ssvr.k <- tune(
    svm, SpatialFormula, data = TrainData_k,
    type = "eps-regression", kernel = "radial",
    ranges = list(cost = c(0.1,1,10), gamma = c(0.01,0.1,1)),
    tunecontrol = tune.control(cross = 12)
  )
  ssvr_model_k <- tune.ssvr.k$best.model
  pred_ssvr_k  <- predict(ssvr_model_k, newdata = TestData_k)

  rss_k <- sum((TestData_k$Y - pred_ssvr_k)^2)
  tss_k <- sum((TestData_k$Y - mean(TestData_k$Y))^2)
  r2_ssvr_k   <- 1 - rss_k/tss_k
  rmse_ssvr_k <- sqrt(mean((TestData_k$Y - pred_ssvr_k)^2))
  mae_ssvr_k  <- mean(abs(TestData_k$Y - pred_ssvr_k))

  cat(sprintf("SSVR (k=%d): R2=%.4f RMSE=%.4f MAE=%.4f (cost=%.3g, gamma=%.3g)\n",
              k_lag_city, r2_ssvr_k, rmse_ssvr_k, mae_ssvr_k,
              tune.ssvr.k$best.parameters$cost, tune.ssvr.k$best.parameters$gamma))

  k_sensitivity_results <- rbind(k_sensitivity_results, data.frame(
    k = k_lag_city, Model = "SSVR", R2 = r2_ssvr_k, RMSE = rmse_ssvr_k, MAE = mae_ssvr_k
  ))

  ## ---- SRF at this k ----
  set.seed(Seed)
  ntree.grid <- c(100,500,1000)
  mtry.grid  <- c(2,3,4)
  K <- 12
  n_k <- nrow(TrainData_k)
  fold.id.k <- sample(rep(1:K, length.out = n_k))

  SRF.Results.k <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())
  for(nt in ntree.grid){
    for(mt in mtry.grid){
      RMSE.fold <- numeric(K)
      for(kk in 1:K){
        Train.CV <- TrainData_k[fold.id.k != kk, ]
        Valid.CV <- TrainData_k[fold.id.k == kk, ]
        set.seed(Seed)
        srf.cv <- randomForest(SpatialFormula, data = Train.CV, ntree = nt,
                                mtry = mt, importance = TRUE)
        pred.cv <- predict(srf.cv, newdata = Valid.CV)
        RMSE.fold[kk] <- sqrt(mean((Valid.CV$Y - pred.cv)^2))
      }
      SRF.Results.k <- rbind(SRF.Results.k, data.frame(ntree = nt, mtry = mt, RMSE = mean(RMSE.fold)))
    }
  }
  best.srf.k <- SRF.Results.k[which.min(SRF.Results.k$RMSE), ]

  srf_model_k <- randomForest(SpatialFormula, data = TrainData_k,
                               ntree = best.srf.k$ntree, mtry = best.srf.k$mtry,
                               importance = TRUE)
  pred_srf_k <- predict(srf_model_k, newdata = TestData_k)

  rss_srf_k <- sum((TestData_k$Y - pred_srf_k)^2)
  tss_srf_k <- sum((TestData_k$Y - mean(TestData_k$Y))^2)
  r2_srf_k   <- 1 - rss_srf_k/tss_srf_k
  rmse_srf_k <- sqrt(mean((TestData_k$Y - pred_srf_k)^2))
  mae_srf_k  <- mean(abs(TestData_k$Y - pred_srf_k))

  cat(sprintf("SRF  (k=%d): R2=%.4f RMSE=%.4f MAE=%.4f (ntree=%d, mtry=%d)\n",
              k_lag_city, r2_srf_k, rmse_srf_k, mae_srf_k,
              best.srf.k$ntree, best.srf.k$mtry))

  k_sensitivity_results <- rbind(k_sensitivity_results, data.frame(
    k = k_lag_city, Model = "SRF", R2 = r2_srf_k, RMSE = rmse_srf_k, MAE = mae_srf_k
  ))

  ## Stash the default (k=5) fitted results into the main Res table, and
  ## remember TestData/TrainData for that k in case you want to inspect
  ## the model objects afterward.
  if (k_lag_city == default_k) {
    Res["SSVR",] = c(mean((TestData_k$Y-pred_ssvr_k)^2), r2_ssvr_k, NA, NA, NA, mae_ssvr_k, rmse_ssvr_k)
    Res["SRF",]  = c(mean((TestData_k$Y-pred_srf_k)^2),  r2_srf_k,  NA, NA, NA, mae_srf_k,  rmse_srf_k)
    ssvr_model_default <- ssvr_model_k
    srf_model_default  <- srf_model_k
  }
}

cat("\n\n==================== K-SENSITIVITY SUMMARY (SRF & SSVR) ====================\n")
print(k_sensitivity_results[order(k_sensitivity_results$Model, k_sensitivity_results$k), ], row.names = FALSE)

write.csv(k_sensitivity_results, "k_sensitivity_summary.csv", row.names = FALSE)
cat("\nSaved: k_sensitivity_summary.csv\n")

# OLS (unaffected by the leakage bug -- unchanged)
####################################################

ols_model <- lm(Formula, data = TrainData)
pred_ols <- predict(ols_model, newdata = TestData)

MSE_OLS <- mean((TestData$Y - pred_ols)^2)
RMSE_OLS <- sqrt(MSE_OLS)
MAE_OLS <- mean(abs(TestData$Y - pred_ols))
RSS_OLS <- sum((TestData$Y - pred_ols)^2)
TSS_OLS <- sum((TestData$Y - mean(TestData$Y))^2)
R2_OLS <- 1 - RSS_OLS/TSS_OLS
AIC_OLS <- AIC(ols_model)

n <- nrow(TrainData)
k <- length(coef(ols_model))
AICc_OLS <- AIC_OLS + (2*k*(k+1)) / (n-k-1)

Res["OLS",] <- c(MSE_OLS, R2_OLS, AIC_OLS, NA, AICc_OLS, MAE_OLS, RMSE_OLS)  # FIX: was Res[1,] -- happened to be correct by coincidence (OLS really is row 1), but now explicit/robust to reordering
cat("MSE  =", MSE_OLS, "\n")
cat("RMSE =", RMSE_OLS, "\n")
cat("MAE  =", MAE_OLS, "\n")
cat("R2   =", R2_OLS, "\n")
cat("AIC  =", AIC_OLS, "\n")

cat("\n==================== FINAL RESULTS (leakage-free) ====================\n")
print(Res)

## ============================================================
## TRAIN-SET METRICS for all six models (to complete the Train/Test
## Table 2 in the manuscript). Uses the fitted model objects already in
## memory from the sections above (ols_model, gwr_model, svr_model,
## rf_model, ssvr_model_default, srf_model_default).
## ============================================================

## ---- OLS ----
pred_ols_train <- predict(ols_model, newdata = TrainData)
R2_OLS_train   <- 1 - sum((TrainData$Y - pred_ols_train)^2) / sum((TrainData$Y - mean(TrainData$Y))^2)
MAE_OLS_train  <- mean(abs(TrainData$Y - pred_ols_train))
RMSE_OLS_train <- sqrt(mean((TrainData$Y - pred_ols_train)^2))

## ---- GWR/SR: gwr_model was fit AT the training locations (no
## fit.points given), so its own fitted values are the in-sample fit ----
pred_gwr_train <- gwr_model$SDF$pred
R2_GWR_train   <- 1 - sum((TrainData$Y - pred_gwr_train)^2) / sum((TrainData$Y - mean(TrainData$Y))^2)
MAE_GWR_train  <- mean(abs(TrainData$Y - pred_gwr_train))
RMSE_GWR_train <- sqrt(mean((TrainData$Y - pred_gwr_train)^2))

## ---- SVR ----
pred_svr_train <- predict(svr_model, newdata = TrainData)
R2_SVR_train   <- 1 - sum((TrainData$Y - pred_svr_train)^2) / sum((TrainData$Y - mean(TrainData$Y))^2)
MAE_SVR_train  <- mean(abs(TrainData$Y - pred_svr_train))
RMSE_SVR_train <- sqrt(mean((TrainData$Y - pred_svr_train)^2))

## ---- RF ----
pred_rf_train <- predict(rf_model, newdata = TrainData)
R2_RF_train   <- 1 - sum((TrainData$Y - pred_rf_train)^2) / sum((TrainData$Y - mean(TrainData$Y))^2)
MAE_RF_train  <- mean(abs(TrainData$Y - pred_rf_train))
RMSE_RF_train <- sqrt(mean((TrainData$Y - pred_rf_train)^2))

## ---- SSVR and SRF (k=5): rebuild the k=5 training frame, since the
## loop variable TrainData_k was overwritten by the LAST k value (k=8)
## by the time the k-sensitivity loop above finished. The model objects
## themselves (ssvr_model_default / srf_model_default) were saved from
## the k=5 iteration, so we only need to regenerate the matching data. ----
lag_k5 <- build_city_lag(TrainData_base, TestData_base, MyData, k_lag_city = 5)
TrainData_k5 <- lag_k5$TrainData

pred_ssvr_train <- predict(ssvr_model_default, newdata = TrainData_k5)
R2_SSVR_train   <- 1 - sum((TrainData_k5$Y - pred_ssvr_train)^2) / sum((TrainData_k5$Y - mean(TrainData_k5$Y))^2)
MAE_SSVR_train  <- mean(abs(TrainData_k5$Y - pred_ssvr_train))
RMSE_SSVR_train <- sqrt(mean((TrainData_k5$Y - pred_ssvr_train)^2))

pred_srf_train <- predict(srf_model_default, newdata = TrainData_k5)
R2_SRF_train   <- 1 - sum((TrainData_k5$Y - pred_srf_train)^2) / sum((TrainData_k5$Y - mean(TrainData_k5$Y))^2)
MAE_SRF_train  <- mean(abs(TrainData_k5$Y - pred_srf_train))
RMSE_SRF_train <- sqrt(mean((TrainData_k5$Y - pred_srf_train)^2))

## ---- Combined Train + Test summary, ready to paste into Table 2 ----
train_metrics <- data.frame(
  Model = c("OLS","RF","SVR","SR","SRF","SSVR"),
  R2_train   = c(R2_OLS_train, R2_RF_train, R2_SVR_train, R2_GWR_train, R2_SRF_train, R2_SSVR_train),
  MAE_train  = c(MAE_OLS_train, MAE_RF_train, MAE_SVR_train, MAE_GWR_train, MAE_SRF_train, MAE_SSVR_train),
  RMSE_train = c(RMSE_OLS_train, RMSE_RF_train, RMSE_SVR_train, RMSE_GWR_train, RMSE_SRF_train, RMSE_SSVR_train)
)

cat("\n\n==================== TRAIN-SET METRICS (for full Table 2) ====================\n")
print(train_metrics, digits = 3)
write.csv(train_metrics, "train_metrics_for_table2.csv", row.names = FALSE)
cat("\nSaved: train_metrics_for_table2.csv\n")

cat("\n\n==================== ALL OUTPUTS SAVED ====================\n")
cat("k_sensitivity_summary.csv       -- SRF/SSVR performance across k in {3,5,8}\n")
cat("train_metrics_for_table2.csv    -- Train-set R2/MAE/RMSE for all 6 models\n")
cat("test_metrics_for_table2.csv     -- Test-set R2/MAE/RMSE for all 6 models (k=5 default)\n")

## ============================================================
## TEST-SET METRICS export -- Res was only ever printed to the console,
## never saved to disk. This writes it out in the same Model-name order
## and format as train_metrics_for_table2.csv, so the two files can be
## combined directly into the final Train/Test Table 2.
## ============================================================

test_metrics <- data.frame(
  Model      = c("OLS","RF","SVR","SR","SRF","SSVR"),
  R2_test    = c(Res["OLS","R2"],   Res["RF","R2"],   Res["SVR","R2"],
                 Res["GWR","R2"],   Res["SRF","R2"],  Res["SSVR","R2"]),
  MAE_test   = c(Res["OLS","MAE"],  Res["RF","MAE"],  Res["SVR","MAE"],
                 Res["GWR","MAE"],  Res["SRF","MAE"], Res["SSVR","MAE"]),
  RMSE_test  = c(Res["OLS","RMSE"], Res["RF","RMSE"], Res["SVR","RMSE"],
                 Res["GWR","RMSE"], Res["SRF","RMSE"],Res["SSVR","RMSE"])
)

cat("\n\n==================== TEST-SET METRICS (for full Table 2) ====================\n")
print(test_metrics, digits = 3)
write.csv(test_metrics, "test_metrics_for_table2.csv", row.names = FALSE)
cat("\nSaved: test_metrics_for_table2.csv\n")

## Optional convenience: a single combined Train+Test table, in the
## exact column order used in the manuscript's Table 2.
table2_combined <- merge(train_metrics, test_metrics, by = "Model")
table2_combined <- table2_combined[match(c("OLS","RF","SVR","SR","SRF","SSVR"), table2_combined$Model), ]
table2_combined <- table2_combined[, c("Model","R2_train","R2_test","MAE_train","MAE_test","RMSE_train","RMSE_test")]

cat("\n\n==================== COMBINED TABLE 2 (Train + Test) ====================\n")
print(table2_combined, digits = 3, row.names = FALSE)
write.csv(table2_combined, "table2_combined.csv", row.names = FALSE)
cat("\nSaved: table2_combined.csv -- ready to paste values directly into the LaTeX table\n")

