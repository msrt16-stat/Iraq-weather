## ============================================================
## Standalone script: rebuilds ONLY the k=5 SRF model (the one reported
## in Table 2), fits it once, and extracts variable importance --
## without repeating the full k in {3,5,8} sensitivity loop from
## weather_R_FINAL.R. Much faster if you just need the importance figure.
##
## Run this as a fresh script (new session is fine) -- it reloads the
## data and rebuilds everything it needs from scratch.
## ============================================================

PROJECT_ROOT <- Sys.getenv("IRAQ_WEATHER_ROOT", unset = getwd())  # set by run_all.R; see README
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "outputs", "09_rebuild_srf_and_importance")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)

library(openxlsx)
library(spdep)
library(randomForest)
library(RANN)
library(ggplot2)

Seed <- 123
K_LAG_CITY <- 5

path <- file.path(PROJECT_ROOT, "Dalya5x.xlsx")   # no interactive prompt -- place the data file directly in PROJECT_ROOT
stopifnot("Data file not found -- check DATA_FILE name/location" = file.exists(path))
MyData <- read.xlsx(path)
colnames(MyData) <- c("ID","Y","X1","X2","X3","X4","X5","Latitude","Longitude")
cat("Column order check -- confirm this matches your actual variable meanings:\n")
print(colnames(MyData))

## ---- 75/25 split, same seed as the main pipeline ----
set.seed(Seed)
n_total <- nrow(MyData)
idx_train <- sample(1:n_total, size = round(0.75 * n_total))
TrainData_base <- MyData[idx_train, ]
TestData_base  <- MyData[-idx_train, ]

## ---- build_city_lag(): identical to the main pipeline ----
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

lag_k5 <- build_city_lag(TrainData_base, TestData_base, MyData, K_LAG_CITY)
TrainData_k5 <- lag_k5$TrainData
TestData_k5  <- lag_k5$TestData

SpatialFormula <- Y ~ X1 + X2 + X3 + X4 + X5 + WY + WX1 + WX2 + WX3 + WX4 + WX5

## ---- 12-fold CV grid search (same as the main pipeline, ntree x mtry) ----
set.seed(Seed)
ntree.grid <- c(100, 500, 1000)
mtry.grid  <- c(2, 3, 4)
K <- 12
n_k <- nrow(TrainData_k5)
fold.id <- sample(rep(1:K, length.out = n_k))

grid.results <- data.frame(ntree = integer(), mtry = integer(), RMSE = numeric())
for (nt in ntree.grid) {
  for (mt in mtry.grid) {
    rmse.fold <- numeric(K)
    for (kk in 1:K) {
      Train.CV <- TrainData_k5[fold.id != kk, ]
      Valid.CV <- TrainData_k5[fold.id == kk, ]
      set.seed(Seed)
      rf.cv <- randomForest(SpatialFormula, data = Train.CV, ntree = nt, mtry = mt)
      pred.cv <- predict(rf.cv, newdata = Valid.CV)
      rmse.fold[kk] <- sqrt(mean((Valid.CV$Y - pred.cv)^2))
    }
    grid.results <- rbind(grid.results, data.frame(ntree = nt, mtry = mt, RMSE = mean(rmse.fold)))
  }
}
best <- grid.results[which.min(grid.results$RMSE), ]
cat(sprintf("Best hyperparameters: ntree=%d, mtry=%d\n", best$ntree, best$mtry))

srf_model_default <- randomForest(SpatialFormula, data = TrainData_k5,
                                   ntree = best$ntree, mtry = best$mtry, importance = TRUE)

## Quick sanity check against the Table 2 numbers (should be close to
## R2=0.753, RMSE=0.570 -- if wildly different, something upstream differs)
pred_test <- predict(srf_model_default, newdata = TestData_k5)
r2_test   <- 1 - sum((TestData_k5$Y - pred_test)^2) / sum((TestData_k5$Y - mean(TestData_k5$Y))^2)
rmse_test <- sqrt(mean((TestData_k5$Y - pred_test)^2))
cat(sprintf("Sanity check -- Test R2=%.4f RMSE=%.4f (Table 2 reports R2=0.753, RMSE=0.570)\n", r2_test, rmse_test))

## ============================================================
## Variable importance extraction and plot
## ============================================================
imp_raw <- importance(srf_model_default)
print(imp_raw)

imp_df <- data.frame(
  Variable_raw = rownames(imp_raw),
  Importance   = imp_raw[, "IncNodePurity"]
)

label_map <- c(
  X1  = "Temperature", X2  = "Top of Atmosphere", X3  = "Wind Speed",
  X4  = "Precipitation", X5  = "Soil Moisture",
  WY  = "Spatial Lag: Humidity", WX1 = "Spatial Lag: Temperature",
  WX2 = "Spatial Lag: Top of Atmosphere", WX3 = "Spatial Lag: Wind Speed",
  WX4 = "Spatial Lag: Precipitation", WX5 = "Spatial Lag: Soil Moisture"
)
## NOTE: verify this mapping matches your actual column meanings before trusting the plot.

imp_df$Variable <- label_map[imp_df$Variable_raw]
imp_df$Importance_pct <- 100 * imp_df$Importance / sum(imp_df$Importance)
imp_df$Type <- ifelse(grepl("^Spatial Lag", imp_df$Variable), "Spatial lag", "Direct predictor")

ggplot(imp_df, aes(x = reorder(Variable, Importance_pct), y = Importance_pct, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Direct predictor" = "#1f4e78", "Spatial lag" = "#a6c8e0")) +
  labs(x = "Predictor", y = "Mean Decrease in Node Impurity (%)", fill = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 18, family = "serif", face = "bold"),
        legend.position = "bottom")

ggsave("variable_importance_corrected.pdf", width = 8, height = 6)

cat("\nSaved: variable_importance_corrected.pdf\n")
cat("\nFull importance ranking (for the paper's discussion text):\n")
print(imp_df[order(-imp_df$Importance_pct), c("Variable","Importance_pct")], row.names = FALSE)
