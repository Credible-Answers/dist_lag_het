# ============================================================================
# DEMO: Replicating Section 5.3 - Newspapers and Presidential Turnout
# ============================================================================
# This script demonstrates the dist_lag_het package using the application
# from Section 5.3 of de Chaisemartin and D'Haultfoeuille (2025):
# "Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption"
#
# The application studies the effect of daily newspapers on presidential
# election turnout in US counties from 1868 to 1928.
# ============================================================================

# Load the package
library(dist_lag_het)
library(MASS)

# Load the data
data <- read.csv("for_estim_RC_model.csv", header = TRUE)
cat("Data loaded:", nrow(data), "rows,", ncol(data), "columns\n")

# ============================================================================
# Set up the estimation
# ============================================================================

# Number of lags (K=2 means current treatment + 2 lags affect outcome)
K <- 2

# Column names
group_col <- "cnty90"          # County identifier
deltaY_col <- "changeprestout" # Change in presidential turnout
deltaD_col <- "changedailies"  # Change in number of dailies
D_col <- "numdailies"          # Number of daily newspapers

# Year dummies (columns 18-30, which are y4 through y16)
# With K=2, we start at column 16+K=18
X_cols <- names(data)[18:30]

# Get number of counties
list_county <- unique(data[[group_col]])
nb_c <- length(list_county)
cat("Number of counties:", nb_c, "\n\n")

# ============================================================================
# Point Estimation
# ============================================================================

cat("=== POINT ESTIMATION ===\n")
start_time <- Sys.time()

result <- estim_RC_model_unified(
  K = K,
  data = data,
  group_col = group_col,
  deltaY_col = deltaY_col,
  deltaD_col = deltaD_col,
  D_col = D_col,
  X_cols = X_cols,
  model = "base"
)

end_time <- Sys.time()
cat("\nEstimation completed in:", format(end_time - start_time), "\n\n")

cat("Beta coefficients:\n")
print(round(result$B_hat, 4))
cat("\nNumber of observations used:\n")
print(result$Nobs)

# ============================================================================
# Bootstrap Standard Errors (200 replications)
# ============================================================================

cat("\n=== BOOTSTRAP STANDARD ERRORS (B=200) ===\n")
start_time <- Sys.time()

result_boot <- estim_RC_model_unified(
  K = K,
  data = data,
  group_col = group_col,
  deltaY_col = deltaY_col,
  deltaD_col = deltaD_col,
  D_col = D_col,
  X_cols = X_cols,
  model = "base",
  bootstrap = TRUE,
  B = 200,
  verbose = TRUE
)

end_time <- Sys.time()
bootstrap_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("\nBootstrap completed in:", round(bootstrap_time, 1), "seconds\n\n")

# Display results
summary_RC_model(result_boot)

# ============================================================================
# Results Table
# ============================================================================

cat("\n=== RESULTS TABLE ===\n")
results_table <- data.frame(
  Coefficient = c("beta_0 (contemporaneous)",
                  "beta_1 (lag 1)",
                  "beta_2 (lag 2)"),
  Estimate = round(result_boot$B_hat, 4),
  Std_Error = round(result_boot$se, 4),
  t_stat = round(result_boot$t_stat, 2),
  CI_lower = round(result_boot$ci_lower, 4),
  CI_upper = round(result_boot$ci_upper, 4),
  N_obs = result_boot$Nobs
)
print(results_table, row.names = FALSE)

# ============================================================================
# Interpretation
# ============================================================================

cat("\n=== INTERPRETATION ===\n")
cat("beta_0:", round(result_boot$B_hat[1], 4),
    "- Effect of an additional daily newspaper on turnout in the same period\n")
cat("beta_1:", round(result_boot$B_hat[2], 4),
    "- Effect of an additional daily newspaper on turnout one period later\n")
cat("beta_2:", round(result_boot$B_hat[3], 4),
    "- Effect of an additional daily newspaper on turnout two periods later\n")

cat("\nNote: Bootstrap runtime with 200 replications:", round(bootstrap_time, 1), "seconds\n")
