# ===============================================================
# FILE: train_rf_models.R
# PURPOSE: Train Baseline and Two-Step Random Forest models for
#          OHC anomaly reconstruction using EN4 synthetic data.
# ===============================================================

# -----------------------------
# 1. Load Required Libraries
# -----------------------------
library(ncdf4)       # For reading NetCDF ocean profile data
library(ranger)      # Fast implementation of Random Forests
library(lubridate)   # For date handling
library(dplyr)       # Data wrangling
library(purrr)       # Functional programming (map-style functions)
library(glue)        # String interpolation for messages
library(parallel)    # Detect number of available CPU cores
library(tidyr)       # Data reshaping and missing value handling
library(ggplot2)     # Optional: For future diagnostics and plotting

# -----------------------------
# 2. Extract Layer-1 Profiles from NetCDF
# -----------------------------
extract_layer1_profiles <- function(nc_path) {
  nc <- nc_open(nc_path)                 # Open NetCDF file
  on.exit(nc_close(nc), add = TRUE)      # Ensure it closes on function exit
  
  lat  <- ncvar_get(nc, "ts_lat")       # Profile latitude
  lon  <- ncvar_get(nc, "ts_lon")       # Profile longitude
  ymd  <- ncvar_get(nc, "en4_ymd")      # Year, month, day for each profile
  mask <- ncvar_get(nc, "dohc_mask_by_en4_maxdepth")  # Valid profile mask
  dadohc1 <- ncvar_get(nc, "dadohc")[1, ]  # Layer-1 anomaly target
  
  # Filter to valid Layer-1 profiles only
  idx <- which(mask[1, ] == 1)
  year  <- as.integer(ymd[1, idx])
  month <- as.integer(ymd[2, idx])
  day   <- as.integer(ymd[3, idx])
  date  <- make_date(year, month, day)
  
  # Assemble into a cleaned tibble
  tibble(
    latitude  = as.numeric(lat[idx]),
    longitude = as.numeric(lon[idx]),  
    date      = date,
    dadohc    = as.numeric(dadohc1[idx])
  ) %>%
    drop_na(dadohc) %>%
    filter(is.finite(dadohc))
}

# -----------------------------
# 3. Feature Engineering
# -----------------------------
featurize <- function(df) {
  df %>% mutate(
    doy      = yday(date),                             # Day of year (1–365)
    sin1     = sin(2 * pi * doy / 365.25),             # First harmonic seasonal cycle
    cos1     = cos(2 * pi * doy / 365.25),
    decoyear = decimal_date(date),                     # Continuous numeric year
    lat2     = latitude^2,                             # Quadratic spatial terms
    lon2     = longitude^2,
    lat_lon  = latitude * longitude                    # Interaction term
  )
}

# -----------------------------
# 4. Random Train/Test Split
# -----------------------------
split_data <- function(df, train_frac = 0.8) {
  set.seed(42)                     
  n <- nrow(df)
  train_idx <- sample(n, size = floor(train_frac * n))
  test_idx  <- setdiff(seq_len(n), train_idx)
  list(
    train = df[train_idx, ],
    test  = df[test_idx,  ]
  )
}

# -----------------------------
# 5. Train Baseline RF Model
# -----------------------------
train_baseline_rf <- function(train) {
  nthreads <- max(1, detectCores() - 1) 
  
  
  # Grid search: mtry (variables tried per split) and min.node.size
  grid <- expand.grid(
    mtry = 3:6,
    min.node.size = c(5, 10, 20)
  )
  
  
  # Loop through hyperparameter grid and collect OOB errors
  oob_res <- grid %>%
    mutate(oob_mse = pmap_dbl(grid, function(mtry, min.node.size) {
      ranger(
        dadohc ~ latitude + longitude + decoyear + sin1 + cos1 + lat2 + lon2 + lat_lon,
        data = train,
        mtry = mtry,
        min.node.size = min.node.size,
        num.trees = 100,
        oob.error = TRUE,
        num.threads = nthreads,
        seed = 42
      )$prediction.error
    })) %>% arrange(oob_mse) # Sort by MSE
  
  
  # Best hyperparameter combination
  best <- oob_res[1, ]
  
  
  # Train final model with 500 trees using best parameters
  t0 <- Sys.time()
  final_model <- ranger(
    dadohc ~ latitude + longitude + decoyear + sin1 + cos1 + lat2 + lon2 + lat_lon,
    data = train,
    mtry = best$mtry,
    min.node.size = best$min.node.size,
    num.trees = 500,
    importance = "impurity",
    oob.error = TRUE,
    num.threads = nthreads,
    seed = 42
  )
  t1 <- Sys.time()
  
  
  # Return model and diagnostics
  list(
    model = final_model,
    best_params = list(mtry = best$mtry, min.node.size = best$min.node.size),
    oob_rmse = sqrt(final_model$prediction.error),
    tune_table = oob_res %>% mutate(oob_rmse = sqrt(oob_mse)),
    train_time_s = as.numeric(difftime(t1, t0, units = "secs"))
  )
}

# -----------------------------
# 6. Train Two-Step RF Model
# -----------------------------
train_two_step_rf <- function(train) {
  nthreads <- max(1, detectCores() - 1) # Use available CPU cores
  set.seed(42)
  
  
  # ---- STEP 1: Fit seasonal model using harmonic + spatial terms ----
  grid1 <- expand.grid(mtry = 3:5, min.node.size = c(5, 10, 20))
  oob1 <- grid1 %>%
    mutate(mse = pmap_dbl(grid1, function(mtry, min.node.size) {
      ranger(
        dadohc ~ sin1 + cos1 + latitude + longitude + lat2 + lon2 + lat_lon,
        data = train,
        mtry = mtry,
        min.node.size = min.node.size,
        num.trees = 100,
        oob.error = TRUE,
        num.threads = nthreads,
        seed = 42
      )$prediction.error
    })) %>% arrange(mse)
  
  
  best1 <- oob1[1, ]
  rf1 <- ranger(
    dadohc ~ sin1 + cos1 + latitude + longitude + lat2 + lon2 + lat_lon,
    data = train,
    mtry = best1$mtry,
    min.node.size = best1$min.node.size,
    num.trees = 500,
    importance = "impurity",
    oob.error = TRUE,
    num.threads = nthreads,
    seed = 42
  )
  
  
  # Generate residuals after removing seasonal fit
  train$seasonal_pred <- predict(rf1, data = train)$predictions
  train$residual <- train$dadohc - train$seasonal_pred
  
  
  # ---- STEP 2: Fit residual model using remaining predictors ----
  grid2 <- expand.grid(mtry = 3:5, min.node.size = c(5, 10, 20))
  oob2 <- grid2 %>%
    mutate(mse = pmap_dbl(grid2, function(mtry, min.node.size) {
      ranger(
        residual ~ decoyear + latitude + longitude + lat2 + lon2 + lat_lon,
        data = train,
        mtry = mtry,
        min.node.size = min.node.size,
        num.trees = 100,
        oob.error = TRUE,
        num.threads = nthreads,
        seed = 42
      )$prediction.error
    })) %>% arrange(mse)
  
  
  best2 <- oob2[1, ]
  rf2 <- ranger(
    residual ~ decoyear + latitude + longitude + lat2 + lon2 + lat_lon,
    data = train,
    mtry = best2$mtry,
    min.node.size = best2$min.node.size,
    num.trees = 500,
    importance = "impurity",
    oob.error = TRUE,
    num.threads = nthreads,
    seed = 42
  )
  
  
  list(
    seasonal_model = rf1,
    residual_model = rf2,
    tune1 = oob1 %>% mutate(oob_rmse = sqrt(mse)),
    tune2 = oob2 %>% mutate(oob_rmse = sqrt(mse)),
    oob_rmse_seasonal = sqrt(rf1$prediction.error),
    oob_rmse_residual = sqrt(rf2$prediction.error)
  )
}


# -----------------------------
# 7. Diagnostic Functions
# -----------------------------

# -- Baseline Model Diagnostics --
print_baseline_diagnostics <- function(test, model) {
  pred <- predict(model, data = test)$predictions
  e <- pred - test$dadohc
  
  rmse <- sqrt(mean(e^2, na.rm = TRUE)) # RMSE
  mae  <- mean(abs(e), na.rm = TRUE) # MAE
  bias <- mean(e, na.rm = TRUE) # Bias
  r2   <- suppressWarnings(cor(pred, test$dadohc, use = "complete.obs")^2) # R-squared
  
  oob_rmse <- if (!is.null(model$prediction.error)) sqrt(model$prediction.error) else NA_real_
  
  cat(glue(
    "Baseline RF - Test RMSE: {round(rmse, 3)}, R²: {round(r2, 3)}\n",
    "OOB RMSE: {round(oob_rmse, 3)}, MAE: {round(mae, 3)}, Bias: {round(bias, 3)}\n"
  ))
  
  return(list(
    test_aug = mutate(test, pred = pred, error = e),
    metrics  = list(rmse = rmse, mae = mae, bias = bias, r2 = r2, oob_rmse = oob_rmse)
  ))
}


# -- Two-Step Model Diagnostics --
print_two_step_diagnostics <- function(test, two_step) {
  test$seasonal_pred <- predict(two_step$seasonal_model, data = test)$predictions
  test$residual_pred <- predict(two_step$residual_model, data = test)$predictions
  test$combined_pred <- test$seasonal_pred + test$residual_pred
  
  rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
  mae  <- function(e) mean(abs(e), na.rm = TRUE)
  bias <- function(e) mean(e, na.rm = TRUE)
  safe_cor2 <- function(x, y) {
    r <- suppressWarnings(cor(x, y, use = "complete.obs"))
    if (is.na(r)) NA_real_ else r^2
  }
  
  e_comb <- test$combined_pred - test$dadohc
  e_seas <- test$seasonal_pred - test$dadohc
  e_res_only <- test$residual_pred - (test$dadohc - test$seasonal_pred)
  
  metrics <- list(
    rmse_combined = rmse(e_comb),
    mae_combined  = mae(e_comb),
    bias_combined = bias(e_comb),
    r2_combined   = safe_cor2(test$combined_pred, test$dadohc),
    
    rmse_seasonal = rmse(e_seas),
    r2_seasonal   = safe_cor2(test$seasonal_pred, test$dadohc),
    
    rmse_residual_only = rmse(e_res_only),
    r2_residual_only   = safe_cor2(test$residual_pred, test$dadohc - test$seasonal_pred)
  )
  
  cat(glue(
    "Two-Step RF - Test RMSE: {round(metrics$rmse_combined, 3)}, ",
    "R²: {round(metrics$r2_combined, 3)}\n"
  ))
  
  return(list(test_aug = test, metrics = metrics))
}


# -----------------------------
# 8. Main Training and Evaluation Workflow
# -----------------------------

# List all NetCDF files in the working directory
nc_files <- list.files(pattern = "\\.nc$")

cat("Reading NetCDF files...\n")
# Read and extract Layer-1 profiles from all NetCDF files
df_raw <- map_dfr(nc_files, extract_layer1_profiles)

cat("Engineering features...\n")
# Apply feature engineering: seasonal harmonics, decimal year, quadratic terms, etc.
df_feat <- featurize(df_raw)

cat("Splitting data into training and test sets...\n")
# Randomly split into training (80%) and test (20%) sets
splits <- split_data(df_feat)

# -----------------------------
# Model Training
# -----------------------------

cat("Training Baseline Random Forest...\n")
baseline <- train_baseline_rf(splits$train)

cat("Training Two-Step Random Forest (Seasonal + Residual)...\n")
two_step <- train_two_step_rf(splits$train)

# -----------------------------
# Model Evaluation (Diagnostics)
# -----------------------------

cat("Evaluating Baseline RF on test set...\n")
diag_baseline <- print_baseline_diagnostics(splits$test, baseline$model)

cat("Evaluating Two-Step RF on test set...\n")
diag_two_step <- print_two_step_diagnostics(splits$test, two_step)

# -----------------------------
# Summary Tables for RMSE, R², OOB error
# -----------------------------

# Extract OOB RMSE from Baseline model
oob_baseline <- if (!is.null(baseline$model$prediction.error)) {
  sqrt(baseline$model$prediction.error)
} else {
  NA_real_
}

# Build main comparison table for Baseline and Two-Step combined prediction
main_tbl <- data.frame(
  model     = c("Baseline RF", "Two-Step RF (combined)"),
  rmse_test = c(diag_baseline$metrics$rmse,
                diag_two_step$metrics$rmse_combined),
  r2_test   = c(diag_baseline$metrics$r2,
                diag_two_step$metrics$r2_combined),
  oob_rmse  = c(oob_baseline, NA_real_)  # Combined model has no single OOB RMSE
)

print(main_tbl)

# Extract component OOB RMSEs
oob_seasonal <- sqrt(two_step$seasonal_model$prediction.error)
oob_residual <- sqrt(two_step$residual_model$prediction.error)

# Build diagnostics table for individual stages of Two-Step RF
component_tbl <- data.frame(
  stage       = c("Seasonal-only model", "Residual-only model"),
  rmse_stage  = c(diag_two_step$metrics$rmse_seasonal,
                  diag_two_step$metrics$rmse_residual_only),
  r2_stage    = c(diag_two_step$metrics$r2_seasonal,
                  diag_two_step$metrics$r2_residual_only),
  oob_rmse    = c(oob_seasonal, oob_residual)
)

print(component_tbl)

# -----------------------------
# Variable Importance
# -----------------------------
cat("Baseline model variable importance:\n")
print(baseline$model$variable.importance)

cat("Two-Step model (seasonal) variable importance:\n")
print(two_step$seasonal_model$variable.importance)

cat("Two-Step model (residual) variable importance:\n")
print(two_step$residual_model$variable.importance)

# -----------------------------
# Save all models and diagnostics
# -----------------------------
saveRDS(baseline,        "baseline_model.rds")
saveRDS(two_step,        "two_step_model.rds")
saveRDS(splits,          "splits_data.rds")
saveRDS(diag_baseline,   "diag_baseline.rds")
saveRDS(diag_two_step,   "diag_two_step.rds")

