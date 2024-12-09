rm(list = ls()) 


load("citibike.RData")

library(dplyr)
library(lubridate)
cb_df <- as.data.frame(citibike)%>%
  mutate(datetime = make_datetime(year, month, day, hour))

library(tsibble)
cb_ts <- as_tsibble(cb_df)

# Seasonality:
library(forecast)
findfrequency(cb_ts) # Frequency: 24 -> Daily seasonality


# Plots
library(ggplot2)
library(papaja)
ggplot(cb_df, aes(x = datetime)) +
  geom_smooth(aes(y = demand, color = "demand")) +
  labs(title = "Demand for Bikes (Period January to May 2023)",
       x = "Date",
       y = "Demand") + 
  scale_color_manual(values = c("demand" = "seagreen"))+
  theme_apa()

ggplot(cb_df, aes(x = as.POSIXct(paste(year, month, day, hour, sep = "-"), format = "%Y-%m-%d-%H"), y = demand)) +
  geom_line(color = "blue") +
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month") + # Adds month names (Jan, Feb, etc.)
  labs(title = "Hourly Bike Demand Over Time", x = "Time (Months)", y = "Demand")+
  theme_apa()


## Demand:weekly 
ggplot(cb_df %>% group_by(wkday) %>% summarize(avg_demand = mean(demand, na.rm = TRUE)), 
       aes(x = wkday, y = avg_demand)) +
  geom_line(color = "purple", size = 1.15) +
  scale_x_continuous(breaks = 1:7, labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  labs(title = "Average Demand by Day of the Week", x = "Day of Week", y = "Average Demand") +
  theme_apa()

### Demand: hourly
ggplot(cb_df, aes(x = factor(hour), y = demand)) +
  geom_boxplot() +
  labs(title = "Hourly Bike Demand", x = "Hour of the Day", y = "Demand" ) +
  theme_apa()


###--------------------------------------------------------------------------###
## Training data (January - April) and test data (may)
trainingdata <- cb_ts %>%
  filter(month <= 4) 

testdata <- cb_ts %>%
  filter(month == 5) 

###---------------------------------------------------------------------------###
# ARIMA:
ggtsdisplay(trainingdata$demand, lag.max = 48)

# Seasonally differenced 
library(feasts)
ggtsdisplay(difference(trainingdata$demand, lag = 24), lag = 48)

# auto select model and recommend model based on ACF/PACF
library(fable)
fit <- trainingdata %>%
  model(ARIMA(demand ~ pdq(2, 0, 0) + PDQ(2, 1, 0)), 
        auto = ARIMA(demand, stepwise = FALSE, approx = FALSE))

library(tidyr)
fit |> pivot_longer(everything(), names_to = "Model name",
                    values_to = "Orders")
# Automatic search suggests ARIMA(1,0,0)(2,1,2)[24]

glance(fit) |> arrange(AICc) 


SARIMA1 <- Arima(trainingdata$demand, order = c(1, 0, 0), seasonal = list(order = c(2, 1, 2), period = 24))
SARIMA2 <- Arima(trainingdata$demand, order = c(2, 0, 0), seasonal = list(order = c(2, 1, 0), period = 24))

      
# Out-of-sample Forecast Exercise

## Expanding window forecast. 
N <- length(cb_ts$demand)
Nstart <- length(trainingdata$demand)
fchorizon <- 24
forecast <- numeric(N)  
forecast[1:Nstart] <- trainingdata$demand

for (i in seq(Nstart, N - fchorizon, by = fchorizon)) {
  
  training <- cb_ts$demand[1:i]  
  
  model <- Arima(training, order = c(1, 0, 0), seasonal = list(order = c(2, 1, 2), period = 24))
  
  forecast[(i + 1):(i + fchorizon)] <- forecast(model, h = fchorizon)$mean
}

testdata$forecast <- forecast[(Nstart+1):N]

ggplot(testdata, aes(x = datetime)) +
  geom_line(aes(y = demand, color = "real demand")) +
  geom_line(aes(y = forecast, color = "expw forecast")) +
  labs(title = "Demand for Bikes (May 2023): SARIMA(1,0,0)(2,1,2)24",
       x = "Date",
       y = "Demand") + 
  scale_color_manual(values = c("real demand" = "seagreen", "expw forecast" = "chocolate1")) +
  theme_apa() 

## Metrics
expw_metrics <- data.frame(RMSE = sqrt(mean((testdata$forecast - testdata$demand)^2)),
                           MAE = mean(abs(testdata$forecast - testdata$demand)),
                           MAPE = mean(abs((testdata$forecast - testdata$demand) / testdata$demand)) * 100)


## Rolling Window forecast. 
rw_forecast <- numeric(nrow(testdata))  


for (i in seq(1, nrow(testdata), by = fchorizon)) {
  
  training <- cb_ts$demand[i:(i + nrow(trainingdata) - 1)] 
  
  rw_model <- Arima(training, order = c(1, 0, 0), seasonal = list(order = c(2, 1, 2), period = 24))
  
  rw_forecast[i:(i+23)] <- forecast(rw_model, h = fchorizon)$mean
  
}

testdata$rw_forecast <- rw_forecast

ggplot(testdata, aes(x = datetime)) +
  geom_line(aes(y = demand, color = "real demand")) +
  geom_line(aes(y = rw_forecast, color = "rw forecast")) +
  labs(title = "Demand for Bikes (May 2023): SARIMA(1,0,0)(2,1,2)24",
       x = "Date",
       y = "Demand") + 
  scale_color_manual(values = c("real demand" = "seagreen", "rw forecast" = "chocolate1")) +
  theme_apa() 

# Metrics
rw_metrics <- data.frame(RMSE = sqrt(mean((testdata$rw_forecast - testdata$demand)^2)),
                         MAE = mean(abs(testdata$rw_forecast - testdata$demand)),
                         MAPE = mean(abs((testdata$rw_forecast - testdata$demand) / testdata$demand)) * 100)

###---------------------------------------------------------------------------###
# EXPONENTIAL SMOOTHING MODELS 

#Estimate appropriate models 

# STL Decomposition
citibike_ts <- ts(trainingdata$demand, start = c(2023, 1), frequency = 24)
decomp <- stl(citibike_ts, s.window = "periodic")
autoplot(decomp, main = "STL Decomposition: Trend, Seasonality, and Residuals")

# Turn the demand into time series
train <- trainingdata %>% select(demand, datetime)

# Fit the model (M, Ad, M)
# damped trend and multiplicative seasonality
library(fpp3)
fit <- train |> 
  model(ANN = ETS(demand ~ error("M")+ trend("Ad") + season ("M")))  


fit2<- train |> 
  model(ANN = ETS(demand))  #I tried the auto-model as well, give A, N, A model

# Print estimated parameters
tidy(fit) 
report(fit)    # model (M, damped, M)
report(fit2)   #model(A, N, A), lower AIC so maybe better? 

#Smoothing parameters for (M, damped, M:
#  alpha = 0.7177319 -> prioritizes recent observations but there is still significant weight assign to past data
#beta  = 0.0002529429 -> stable, weak trend (changes very gradually)
#gamma = 0.01643204 -> stable seasonality, seasonal patterns are consistent and do not change dramatically over time
#phi   = 0.9799506 -> minimal trend damping, continue current direction strongly in future

# Forecast
fc <- fit |> 
  forecast(h = 48)  #originally 24 but i tried with 48 to see the seasonality

fc2 <- fit2 |>
  forecast(h=48)


# Plot forecast (full train + forecast)
fc |>    #fit1
  autoplot(train) +
  geom_line(aes(y = .fitted), col = "#D55E00", data = augment(fit)) +
  labs(
    y = "Demand Forecast",
    x = "Datetime",
    title = "Demand Forecast Using Holt-Winters Method (M, damped, M)"
  ) +
  guides(colour = "none")

#Plot only the forecast
fc |>
  autoplot() +  
  labs(
    y = "Demand Forecast",
    x = "Datetime",
    title = "Demand Forecast Using Holt-Winters Method (M, Damped, M)"
  ) +
  guides(colour = "none")


# Plot fitted train data and forecast with model 2
fc2 |>   #fit2
  autoplot(train) +
  geom_line(aes(y = .fitted), col = "#D55E00", data = augment(fit2)) +
  labs(
    y = "Demand Forecast",
    x = "Datetime",
    title = "Demand Forecast Using Holt-Winters Method (A, N, A)"
  ) +
  guides(colour = "none")

#Plot only the forecast
fc2 |>
  autoplot() +  # Plot only the forecast
  labs(
    y = "Demand Forecast",
    x = "Datetime",
    title = "Demand Forecast Using Holt-Winters Method (A, N, A)"
  ) +
  guides(colour = "none")

  
#Print the forecast 
fc_table <- fc |> as_tibble()
print(fc_table, n = Inf, width = Inf)  

fc_table2 <- fc2 |> as_tibble()
print(fc_table2, n = Inf, width = Inf)  
# The demand forecasted is a distribution N(mean, variance)



# Exponential Smoothing: Out-of-sample forecast

# Function to calculate error metrics
calculate_metrics <- function(actual, predicted, model_name) {
  tibble(
    Model = model_name,
    RMSE = sqrt(mean((predicted - actual)^2, na.rm = TRUE)),
    MAE = mean(abs(predicted - actual), na.rm = TRUE),
    MAPE = mean(abs((predicted - actual) / actual), na.rm = TRUE) * 100
  )
}

# Initialize storage for forecasts and metrics
store_forecasts <- function(forecast, validation_window, model_name, storage_list) {
  filtered_forecast <- forecast %>% filter(datetime %in% validation_window$datetime)
  storage_list[[model_name]] <- append(storage_list[[model_name]], list(filtered_forecast))
  storage_list
}



# Rolling Window
train_rows <- nrow(train)
forecast_horizon <- 24
num_windows <- nrow(cb_ts) - train_rows

rolling_forecasts <- list(MADM = list(), ANA = list())
rolling_metrics <- list()

for (i in seq(1, num_windows, by = forecast_horizon)) {
  train_window <- cb_ts[i:(i + train_rows - 1), ]
  validation_window <- cb_ts[(i + train_rows):(i + train_rows + forecast_horizon - 1), ]
  
  fit_madm <- train_window %>% model(ETS_MAM = ETS(demand ~ error("M") + trend("Ad") + season("M")))
  fit_ana <- train_window %>% model(ETS_ANA = ETS(demand ~ error("A") + trend("N") + season("A")))
  
  forecast_madm <- fit_madm %>% forecast(h = forecast_horizon)
  forecast_ana <- fit_ana %>% forecast(h = forecast_horizon)
  
  actual <- validation_window$demand
  predicted_madm <- forecast_madm %>% filter(datetime %in% validation_window$datetime) %>% pull(.mean)
  predicted_ana <- forecast_ana %>% filter(datetime %in% validation_window$datetime) %>% pull(.mean)
  
  rolling_metrics[[i]] <- bind_rows(
    calculate_metrics(actual, predicted_madm, "MADM"),
    calculate_metrics(actual, predicted_ana, "ANA")
  )
  
  rolling_forecasts <- store_forecasts(forecast_madm, validation_window, "MADM", rolling_forecasts)
  rolling_forecasts <- store_forecasts(forecast_ana, validation_window, "ANA", rolling_forecasts)
}




# Expanding Window
expanding_forecasts <- list(MADM = list(), ANA = list())
expanding_metrics <- list()

for (i in seq(1, num_windows, by = forecast_horizon)) {
  train_window <- cb_ts[1:(train_rows + i - 1), ]
  validation_window <- cb_ts[(train_rows + i):(train_rows + i + forecast_horizon - 1), ]
  
  fit_madm <- train_window %>% model(ETS_MAM = ETS(demand ~ error("M") + trend("Ad") + season("M")))
  fit_ana <- train_window %>% model(ETS_ANA = ETS(demand ~ error("A") + trend("N") + season("A")))
  
  forecast_madm <- fit_madm %>% forecast(h = forecast_horizon)
  forecast_ana <- fit_ana %>% forecast(h = forecast_horizon)
  
  actual <- validation_window$demand
  predicted_madm <- forecast_madm %>% filter(datetime %in% validation_window$datetime) %>% pull(.mean)
  predicted_ana <- forecast_ana %>% filter(datetime %in% validation_window$datetime) %>% pull(.mean)
  
  expanding_metrics[[i]] <- bind_rows(
    calculate_metrics(actual, predicted_madm, "MADM"),
    calculate_metrics(actual, predicted_ana, "ANA")
  )
  
  expanding_forecasts <- store_forecasts(forecast_madm, validation_window, "MADM", expanding_forecasts)
  expanding_forecasts <- store_forecasts(forecast_ana, validation_window, "ANA", expanding_forecasts)
}

# Combine Metrics
rolling_metrics_df <- bind_rows(rolling_metrics)
expanding_metrics_df <- bind_rows(expanding_metrics)

# Summary Metrics
rolling_summary <- rolling_metrics_df %>%
  group_by(Model) %>%
  summarise(across(RMSE:MAPE, \(x) mean(x, na.rm = TRUE)))
expanding_summary <- expanding_metrics_df %>% group_by(Model) %>% summarise(across(RMSE:MAPE, mean, na.rm = TRUE))

# Print Summary Metrics
print(rolling_summary)
print(expanding_summary)

# Plots
plot_forecasts <- function(forecast_df, actual_df, title, color) {
  ggplot() +
    geom_line(data = forecast_df, aes(x = datetime, y = .mean), color = color, size = 1) +
    geom_line(data = actual_df, aes(x = datetime, y = demand), color = "black", linetype = "dashed") +
    labs(title = title, x = "Datetime", y = "Demand") +
    theme_apa()
}

rolling_forecast_datetimes <- bind_rows(lapply(rolling_forecasts$MADM, as_tibble)) %>%
  pull(datetime)
actual_rolling <- cb_ts %>%
  filter(datetime %in% rolling_forecast_datetimes)
actual_rolling <- as_tsibble(actual_rolling, index = datetime)

expanding_forecast_datetimes <- bind_rows(lapply(expanding_forecasts$MADM, as_tibble)) %>%
  pull(datetime)
actual_expanding <- cb_ts %>%
  filter(datetime %in% expanding_forecast_datetimes)
actual_expanding <- as_tsibble(actual_expanding, index = datetime)

plot_forecasts(bind_rows(rolling_forecasts$MADM), actual_rolling, "Rolling Window - M, Ad,M", "blue")
plot_forecasts(bind_rows(rolling_forecasts$ANA), actual_rolling, "Rolling Window - A, N, A", "red")
plot_forecasts(bind_rows(expanding_forecasts$MADM), actual_expanding, "Expanding Window - M, Ad, M", "blue")
plot_forecasts(bind_rows(expanding_forecasts$ANA), actual_expanding, "Expanding Window - A, N, A", "red")

# Extract ANA forecasts and calculate errors for rolling and expanding windows
errors_ets_rw <- lapply(rolling_forecasts$ANA, function(forecast) {
  forecast_df <- as_tibble(forecast)
  actual_df <- cb_ts %>% filter(datetime %in% forecast_df$datetime) %>% pull(demand)
  actual_df - forecast_df$.mean
})

errors_ets_expw <- lapply(expanding_forecasts$ANA, function(forecast) {
  forecast_df <- as_tibble(forecast)
  actual_df <- cb_ts %>% filter(datetime %in% forecast_df$datetime) %>% pull(demand)
  actual_df - forecast_df$.mean
})

###---------------------------------------------------------------------------###
# TESTS

# Extracted ETS errors for Rolling and Expanding Windows
errors_ets_rw <- unlist(errors_ets_rw)  # Flattened errors from ETS rolling forecasts
errors_ets_expw <- unlist(errors_ets_expw)  # Flattened errors from ETS expanding forecasts

# Replace SARIMA error extraction with actual SARIMA forecasts
errors_sarima_rw <-  testdata$demand - testdata$rw_forecast # Rolling Window SARIMA errors
errors_sarima_expw <- testdata$demand - testdata$forecast  # Expanding Window SARIMA errors


# Diebold-Mariano Test for Rolling Window
dm_rw_result <- dm.test(
  errors_ets_rw,         
  errors_sarima_rw,      
  alternative = "two.sided",
  h = 24,                   
  power = 2,
  varestimator = "acf"
)

# Diebold-Mariano Test for Expanding Window
dm_expw_result <- dm.test(
  errors_ets_expw,       
  errors_sarima_expw,    
  alternative = "two.sided",
  h = 24,                   
  power = 2,
  varestimator = "acf"
)

# Model Confidence Set for Rolling Window
library(MCS)
errors_matrix_rw <- cbind(
  ETS_RW = errors_ets_rw^2,
  SARIMA_RW = errors_sarima_rw^2
)

mcs_rw_result <- MCSprocedure(
  Loss = errors_matrix_rw,
  alpha = 0.05,
  B = 1000,
  statistic = "Tmax"
)


# Model Confidence Set for Expanding Window
errors_matrix_expw <- cbind(
  ETS_EXPW = errors_ets_expw^2,
  SARIMA_EXPW = errors_sarima_expw^2
)

mcs_expw_result <- MCSprocedure(
  Loss = errors_matrix_expw,
  alpha = 0.05,
  B = 1000,
  statistic = "Tmax"
)


