rm(list = ls())
if (!require("forecast")) install.packages("forecast", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
library(forecast)
library(ggplot2)

# Set working directory and read data
setwd("C:/Users/YourUsername/Desktop")
df1 <- read.csv("ARIMA-GM(1,1) analysis data.csv", header = TRUE)

# Check if the data contains missing values
if (any(is.na(df1))) {
  stop("There are missing values in the data, please handle them before continuing.")
}

# Extract years and drug resistance rate data
years <- df1[, 1]  # The first column is the year
drug_resistance_rate <- df1[, 7]  # The seventh column is the drug resistance rate data for that year

# GM(1,1) model function
gm1_model <- function(data) {
  n <- length(data)
  if (n < 2) {
    stop("The data length must be greater than 1.")
  }
  
  # Cumulative sum sequence (AGO)
  x1 <- cumsum(data)
  
  # Calculate adjacent mean generating sequence
  x0_ave <- (x1[1:(n - 1)] + x1[2:n]) / 2
  B_matrix <- cbind(-x0_ave, rep(1, n - 1))
  Y <- data[2:n]
  
  # Least squares estimation of a and b
  coefs <- solve(t(B_matrix) %*% B_matrix) %*% t(B_matrix) %*% Y
  a <- coefs[1]
  b <- coefs[2]
  
  # Get predicted values of the time response function
  future_series <- numeric(n)
  future_series[1] <- x1[1]  # The first value remains unchanged
  
  for (i in 2:n) {
    future_series[i] <- (x1[1] - b / a) * exp(-a * (i - 1)) + b / a
  }
  
  # Calculate fitted values
  fitted_values <- future_series
  
  # Calculate residuals and related metrics
  residuals <- data - fitted_values
  avg_real <- mean(data)
  s1_squared <- var(data)  # Variance of the original data
  s2_squared <- var(residuals)  # Variance of the residuals
  
  # Calculate average relative residual
  epsilon_r <- mean(abs(residuals[2:n]) / data[2:n]) * 100
  
  # Calculate posterior variance ratio C; smaller is better
  C <- s2_squared / s1_squared
  
  # Calculate small probability error P
  P <- mean(residuals < 0.6745 * sd(residuals))
  
  return(list(future_series = future_series, a = a, b = b, C = C, P = P, epsilon_r = epsilon_r))
}

# Use GM(1,1) model to process data
gm1_results <- gm1_model(drug_resistance_rate)
accumulated_data <- gm1_results$future_series
a_coefficient <- gm1_results$a
b_coefficient <- gm1_results$b
C_value <- gm1_results$C
P_value <- gm1_results$P
epsilon_r <- gm1_results$epsilon_r

# Print the calculated results
cat("\nGM(1,1) model estimation results:\n")
cat("Estimated a:", round(a_coefficient, 3), "\n")
cat("Estimated b:", round(b_coefficient, 3), "\n")
cat("Small probability error P:", round(P_value, 3), "\n")

# Fit using auto.arima
garma_model <- auto.arima(drug_resistance_rate)

# Future value prediction
forecasted_values <- forecast(garma_model, h = 5)
future_years <- (max(years) + 1):(max(years) + 5)

# Create data frame
forecasted_data <- data.frame(
  Year = c(years, future_years), 
  Values = c(drug_resistance_rate, forecasted_values$mean)
)

# Create confidence interval data frame
confidence_intervals <- data.frame(
  Year = future_years,
  Lower = forecasted_values$lower[, 2],  # 95% lower bound
  Upper = forecasted_values$upper[, 2],  # 95% upper bound
  Mean = forecasted_values$mean  # Predicted mean
)

# Residual checks
checkresiduals(garma_model)

# Calculate and print RMSE of ARIMA model
arima_predictions <- forecast(garma_model, h = length(drug_resistance_rate))$mean  # Prediction for known data points
arima_rmse <- sqrt(mean((drug_resistance_rate - arima_predictions)^2))  # Calculate RMSE for ARIMA prediction
cat("ARIMA model root mean square error (RMSE):", round(arima_rmse, 3), "\n")

# ACF and PACF plots
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # Set layout and margins
Acf(residuals(garma_model), main = "ACF of ARIMA Residuals")
Pacf(residuals(garma_model), main = "PACF of ARIMA Residuals")

# Print residuals and conduct ACF, PACF analysis
residuals <- residuals(garma_model)

# Set multipanel layout
par(mfrow = c(3, 1), mar = c(5, 5, 3, 1))  # 3 rows 1 column, set margins

# Plot standardized residuals
plot(residuals, type = "h", col = "black", main = "Standardized Residual", xlab = "Time", ylab = "Residual")
abline(h = 0, col = "blue", lwd = 2)  # Add zero line

# Plot ACF
Acf(residuals, main = "ACF of Residuals", lag.max = 20)

# Plot p-values of Ljung-Box statistic
Box.test(residuals, lag = 20, type = "Ljung-Box")  # Ljung-Box test results
ljung_box_stats <- sapply(1:10, function(lag) Box.test(residuals, lag = lag, type = "Ljung-Box")$p.value)

# Plot p-value graph
plot(ljung_box_stats, type = "b", ylim = c(0, 1), xlab = "Lag", ylab = "p value", 
     main = "p value for Ljung-Box Statistic", pch = 19, col = "black")
abline(h = 0.05, col = "blue", lty = 2)  # Add significance level line

# Restore default graphics parameters
par(mfrow = c(1, 1))

# Print AIC and BIC of the optimal ARIMA model
cat("\nOptimal ARIMA model:\n")
cat("AIC:", AIC(garma_model), "\n")
cat("BIC:", BIC(garma_model), "\n")

# Create plot and add labels
forecast_plot <- ggplot() +
  geom_line(aes(x = years, y = drug_resistance_rate, color = "Actual Values"), linewidth = 1) +
  geom_line(aes(x = years, y = fitted(garma_model), color = "Fitted Values"), linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = future_years, y = confidence_intervals$Mean, color = "Forecasted Values"), linewidth = 1, linetype = "dashed") +
  geom_point(aes(x = years, y = drug_resistance_rate, color = "Actual Values"), size = 3) +
  geom_point(aes(x = years, y = fitted(garma_model), color = "Fitted Values"), size = 3) +
  geom_point(aes(x = future_years, y = confidence_intervals$Mean, color = "Forecasted Values"), size = 3) +
  geom_ribbon(data = confidence_intervals, aes(x = Year, ymin = Lower, ymax = Upper), alpha = 0.2, fill = "green") +  # Add confidence interval
  scale_x_continuous(breaks = seq(min(years), max(years) + 5, by = 1)) +
  scale_y_continuous(limits = c(0, max(c(drug_resistance_rate, forecasted_values$mean)) * 1.1), 
                     breaks = pretty(c(0, max(c(drug_resistance_rate, forecasted_values$mean))), n = 5)) +
  labs(title = "MRSA", x = "Year", y = "Resistance Rate (%)") +  
  scale_color_manual(
    values = c("Actual Values" = "blue", "Fitted Values" = "red", "Forecasted Values" = "green"),
    breaks = c("Actual Values", "Fitted Values", "Forecasted Values")  # Specify legend order
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12, family = "Times New Roman"),
    axis.title.x = element_text(size = 12, family = "Microsoft YaHei"),  
    axis.title.y = element_text(size = 12, family = "Microsoft YaHei"),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, family = "Times New Roman"), 
    axis.text.y = element_text(size = 12, family = "Times New Roman"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.ticks = element_line(linewidth = 0.25, color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1.5, "cm"),  # Set legend key (color block) size
    legend.text = element_text(size = 12),  # Set legend text size
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

# Display plot
print(forecast_plot)

# Save the plot
ggsave("MRSA.JPEG", width = 7, height = 4, dpi = 300)
