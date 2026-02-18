# ===============================================================
# FILE: rf_fullfield.R
# PURPOSE:
#   1) Obs-point diagnostics for Baseline and Two-step RF
#   2) Build full-field grid over North Atlantic (0–60N, 280–360E)
#   3) Predict monthly fields 1993–2014 (already anomalies)
#   4) Ocean-mask, save all outputs for plotting
# NOTES:
#   - ME4OH Experiment A uses detrended, de-seasonalised anomalies.
#   - Do NOT subtract a climatology from predictions when plotting.
#   - T0 and monthly climatology here are DIAGNOSTICS only.
# ===============================================================

# -----------------------
# 0. Settings (North Atlantic example)
# -----------------------
basin_lon_min <- 280
basin_lon_max <- 360
basin_lat_min <-   0
basin_lat_max <-  60
grid_res <- 0.5

# -----------------------
# 1. Libraries
# -----------------------
library(ranger)
library(dplyr)
library(ggplot2)
library(lubridate)
library(glue)
library(purrr)
library(rnaturalearth)
library(sf)
sf::sf_use_s2(FALSE)

# -----------------------
# 2. Load models & data splits
# -----------------------
baseline   <- readRDS("baseline_model.rds")           # expects $model
two_step   <- readRDS("two_step_model.rds")           # expects $seasonal_model, $residual_model
splits     <- readRDS("splits_data.rds")              # expects $train, $test
# if you already have a prebuilt masked grid, you can reuse it, but we rebuild below

# -----------------------
# 3. Obs-point diagnostics (on withheld test set)
# -----------------------
cat("\n[Obs-point diagnostics]\n")

# Baseline RF
pred_base_obs <- predict(baseline$model, data = splits$test)$predictions
rmse_base <- sqrt(mean((pred_base_obs - splits$test$dadohc)^2, na.rm = TRUE))
r2_base   <- suppressWarnings(cor(pred_base_obs, splits$test$dadohc, use = "complete.obs")^2)
cat(glue("Baseline RF  → RMSE: {round(rmse_base,3)}, R²: {round(r2_base,3)}\n"))

# Two-step RF (seasonal + residual)
pred_seas_obs <- predict(two_step$seasonal_model, data = splits$test)$predictions
pred_resi_obs <- predict(two_step$residual_model, data = splits$test)$predictions
pred_two_obs  <- pred_seas_obs + pred_resi_obs
rmse_two <- sqrt(mean((pred_two_obs - splits$test$dadohc)^2, na.rm = TRUE))
r2_two   <- suppressWarnings(cor(pred_two_obs, splits$test$dadohc, use = "complete.obs")^2)
cat(glue("Two-step RF   → RMSE: {round(rmse_two,3)}, R²: {round(r2_two,3)}\n"))


# -----------------------
# 4. Full-field prediction
# -----------------------
make_prediction_for_month <- function(year, month) {
  # Build lon×lat grid at one date (0–360 lon domain)
  expand.grid(
    longitude = seq(basin_lon_min, basin_lon_max, by = grid_res),
    latitude  = seq(basin_lat_min, basin_lat_max, by = grid_res),
    date      = as.Date(glue("{year}-{sprintf('%02d',month)}-01"))
  ) |>
    # Create features; extra cols are safe even if Baseline trained without them
    mutate(
      doy      = yday(date),
      sin1     = sin(2*pi*doy/365.25),
      cos1     = cos(2*pi*doy/365.25),
      decoyear = decimal_date(date),
      lat2     = latitude^2,
      lon2     = longitude^2,
      lat_lon  = latitude * longitude
    ) |>
    # Predict both models (outputs are ANOMALIES by definition)
    mutate(
      baseline_pred = predict(baseline$model, data = .)$predictions,
      two_step_pred = predict(two_step$seasonal_model, data = .)$predictions +
        predict(two_step$residual_model, data = .)$predictions
    )
}

# Sanity check: training columns vs features produced here
in_train  <- names(splits$train)
in_future <- names(make_prediction_for_month(2005, 1))
if (length(setdiff(setdiff(in_train, c("dadohc")), in_future)) > 0) {
  warning("Some training features are missing in make_prediction_for_month()!")
}

# -----------------------
# 5. Compute full-field predictions for 1993–2014
# -----------------------
cat("\n[Building full-field predictions 1993–2014]\n")
dates_all <- seq(as.Date("1993-01-01"), as.Date("2014-12-01"), by = "1 month")

all_preds <- map_dfr(dates_all, function(d) {
  y <- year(d); m <- month(d)
  message(" → predicting ", y, "-", sprintf("%02d", m))
  make_prediction_for_month(y, m)
})

cat(glue("Done. Rows = {format(nrow(all_preds), big.mark = ',')}\n"))

# -----------------------
# 6. Ocean mask (0–360 aware)
# -----------------------
build_ocean_mask <- function(df, lon_col = "longitude", lat_col = "latitude") {
  stopifnot(all(c(lon_col, lat_col) %in% names(df)))
  grid_xy <- distinct(df, !!rlang::sym(lon_col), !!rlang::sym(lat_col)) |>
    rename(lon0360 = !!rlang::sym(lon_col), lat = !!rlang::sym(lat_col)) |>
    mutate(lon180 = ifelse(lon0360 > 180, lon0360 - 360, lon0360))
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |> st_geometry()
  pts  <- st_as_sf(grid_xy, coords = c("lon180", "lat"), crs = 4326, remove = FALSE)
  on_land <- lengths(st_intersects(pts, land)) > 0
  transmute(grid_xy, longitude = lon0360, latitude = lat, is_ocean = !on_land)
}

mask_xy <- build_ocean_mask(all_preds)
all_preds_ocean <- inner_join(all_preds, filter(mask_xy, is_ocean),
                              by = c("longitude", "latitude"))

cat("Grid points (unique) before/after mask: ",
    nrow(distinct(all_preds, longitude, latitude)), " / ",
    nrow(distinct(all_preds_ocean, longitude, latitude)), "\n")

# -----------------------
# 7. Save outputs for plotting
# -----------------------
saveRDS(all_preds,            "all_preds_1993_2014.rds")        # full grid (land + ocean)
saveRDS(all_preds_ocean,      "all_preds_northatlantic.rds")    # ocean-only
saveRDS(mask_xy,              "northatlantic_mask_xy.rds")

cat("\n[Saved]: all_preds_1993_2014.rds, all_preds_northatlantic.rds, northatlantic_mask_xy.rds\n")
