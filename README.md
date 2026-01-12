# dist_lag_het

R package to estimate distributed lag regressions with heterogeneous treatment effects.

Based on the methodology described in Sections 4.2 and 4.3 of *"Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption"* by Clement de Chaisemartin and Xavier D'Haultfoeuille.

## Installation

You can install the package directly from GitHub:

```r
# Install from GitHub using devtools
devtools::install_github("chaisemartinpackages/dist_lag_het")

# Or using remotes
remotes::install_github("chaisemartinpackages/dist_lag_het")


```

## Quick Start

### Estimate Models

#### Using the Unified Interface

```r
library(DistLagHet)
# Base model
result <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  weights = weights,
  model = "base"
)

# View results
summary_RC_model(result)

# Estimate with bootstrap standard errors
result_boot <- estim_RC_model_unified(
  K = 2,
  data = data,
  group_col = "group",
  deltaY_col = "delta_Y",
  deltaD_col = "delta_D",
  D_col = "D",
  X_cols = c("X1", "X2"),
  model = "base",
  bootstrap = TRUE,  # Enable bootstrap
  B = 1000           # Number of bootstrap iterations
)

summary_RC_model(result_boot)
```

## Application: Newspapers and Presidential Turnout (Section 5.3)

This example replicates the application from Section 5.3 of the paper, which studies the effect of the number of daily newspapers on presidential election turnout in US counties from 1868 to 1928.

```r
library(dist_lag_het)

# Load the data
data <- read.csv("for_estim_RC_model.csv", header = TRUE)

# Set up the estimation with K=2 lags
K <- 2
group_col <- "cnty90"          # County identifier
deltaY_col <- "changeprestout" # Change in presidential turnout
deltaD_col <- "changedailies"  # Change in number of dailies
D_col <- "numdailies"          # Number of daily newspapers
X_cols <- names(data)[18:30]   # Year dummies (y4 through y16)

# Get number of counties for weights
list_county <- unique(data[[group_col]])
nb_c <- length(list_county)

# Estimate the model
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

# View point estimates
print(result$B_hat)

# Estimate with bootstrap standard errors (200 reps, ~30-60 seconds)
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
  B = 200
)

# Display results with standard errors
summary_RC_model(result_boot)
```

The model estimates the effect of having an additional daily newspaper on presidential turnout, allowing for effects of current treatment (`beta_0`) and up to two lags (`beta_1`, `beta_2`).

### Results

| Coefficient | Estimate | Std. Error | t-stat | 95% CI |
|-------------|----------|------------|--------|--------|
| beta_0 (contemporaneous) | 0.0042 | 0.0022 | 1.92 | [-0.0001, 0.0085] |
| beta_1 (lag 1) | -0.0014 | 0.0025 | -0.57 | [-0.0063, 0.0035] |
| beta_2 (lag 2) | -0.0016 | 0.0021 | -0.79 | [-0.0057, 0.0024] |

*Runtime: ~12 seconds for 200 bootstrap replications*

**Data**: 1,147 counties, 16,249 observations

## Data Format

Your data should be in long format with the following columns:

- **group identifier**: Column identifying groups (e.g., counties, individuals)
- **D**: Treatment level
- **delta_D**: First difference of treatment (ΔD)
- **Y**: Outcome variable
- **delta_Y**: First difference of outcome (ΔY)
- **X1, X2, ...**: Covariates

The first differences should exclude the first period (which would have NA values).

## Options

K: the number of treatment lags assumed to affect the current outcome. For instance, K=2 assumes that the current treatment and its first two lags affect the outcome.

The package provides three model specifications:
1. **Base Model** (`model = "base"`): Standard distributed-lag model with K lags
2. **Interactions** (`model = "interactions"`): Includes interaction terms between current treatments and lags
3. **Full Dynamics** (`model = "full_dynamics"`): Allows all treatment lags up to period one to affect the outcome (then K does not need to be specified).

weights: to weight the estimation.

## Functions

### Main Estimation Functions

- `estim_RC_model_unified()`: Unified interface for all models with optional bootstrap standard errors

### Utility Functions

- `summary_RC_model()`: Print formatted summary of results

## Output Structure

All estimation functions return a list with:

- `gamma`: Coefficients for covariates
- `B_hat`: Treatment effect coefficients (betas)
- `Nobs`: Number of observations used for each coefficient
- `model`: Model type ("base", "full_dynamics", or "interactions")


## Dependencies

- R (>= 3.5.0)
- MASS
- Rcpp (>= 1.0.0)
- RcppEigen (for C++ optimizations)

## License

MIT License - see LICENSE file for details

## Authors

- Clement de Chaisemartin
- Xavier D'Haultfoeuille
- Henri Fabre (package maintainer)

## Citation

If you use this package in your research, please cite:

```
de Chaisemartin, Clement and d'Haultfoeuille, Xavier, Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption (September 04, 2025). Available at SSRN: https://ssrn.com/abstract=5440734 or http://dx.doi.org/10.2139/ssrn.5440734
```

## Support

For issues and questions:
- Open an issue on GitHub
- Contact the package maintainer : henri.fabre@ensae.fr

## References

de Chaisemartin, Clement and d'Haultfoeuille, Xavier, Treatment-Effect Estimation in Complex Designs under a Parallel-trends Assumption (September 04, 2025). Available at SSRN: https://ssrn.com/abstract=5440734 or http://dx.doi.org/10.2139/ssrn.5440734
