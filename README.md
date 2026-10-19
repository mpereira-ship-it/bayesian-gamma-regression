# Bayesian Gamma regression — College application

R code for the empirical application to `ISLR::College` (777 observations, `Outstate` response). The repository contains HB-Lasso, HB-Ridge, SS-BEN, and the three Gamma comparators.

Run from the repository root:

```r
install.packages(c("ISLR", "posterior", "ggplot2", "gamboostLSS", "mboost"))
```

```bash
Rscript run_college.R
```

The frozen training/test design uses seed `20261012` and 17 standardized covariates in both predictors (`tau = 0`). The final runs use four chains per Bayesian method: 60,000 iterations for HB-Lasso and HB-Ridge, and 180,000 for SS-BEN. GammaLSS uses stratified 5-fold cross-validation on `mstop = 20, 40, ..., 2000`.

`results/` contains the reported selection and test-performance tables. Running the script writes fitted chains, diagnostics, predictions, and figures to the ignored `outputs/` directory. Full MCMC reproduction may require substantial computation.

Reference environment: R 4.5.2; ISLR 1.4; posterior 1.7.0; gamboostLSS 2.2-0; mboost 2.9-14; ggplot2 4.0.3.
