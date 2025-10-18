# BGGM 2.1.5.9000 (development version)
### Major changes to ordinal sampler
- **Stan-style latent centering for ordinal models**: Both the Albert and Cowles ordinals samplers have been refactored to improve numerical stability and mixing in the presence of skewed or ceiling/floor ordinal items:
  - Thresholds are now initialized from the empirical category frequencies of each variable (i.e., cumulative proportions mapped to the probit scale), rather than arbitrary or equally spaced cut-points.
  - Latent variable draws \(Z\) are initialized at the expected value of the truncated normal (conditional mean), rather than uniform draws across the truncation interval. This ensures that ceiling or floor items start in the correct region of the latent space.
  - After each update/sweep of a latent column \(Z_j\), the column is recentered (mean ≈ 0) and the corresponding thresholds are shifted by the same offset, preserving the likelihood but aligning the latent origin. This mirrors the identification scheme used in Stan’s ordered-probit/ordered-logit models (latent mean fixed at 0, thresholds floating).
  - Probit‐scale bounds (±8 on the standard normal scale) have been introduced to cap semi-infinite truncation regions (e.g., \((\tau_{K-1},∞)\)) to avoid numerical overflow and improve stability for extreme category distributions.

### Bug fixes & improvements
- Synchronized threshold matrices in the Cowles sampler: the `current_thresh`, `candidate_thresh`, `thresh_mat`, `c_thresh_mat`, and `thresh_mcmc` now begin from the same baseline after initialization to avoid drift in Metropolis proposals.
- Improved sampler performance for ordinal variables with heavy tails or heavy ceiling/floor effects. This should lead to reduced shrinkage of partial-correlations toward zero when items are highly skewed.

### Compatibility notes
- The statistical model remains unchanged: you still get the same posterior for latent precision/correlation matrices. The changes are purely in parameterisation and initialization of the latent \(Z\) and threshold variables.

# BGGM 2.1.4
## Bug Fixes and Improvements

### C++ `search` Function
1. **Improved Initialization**  
   - The initial adjacency matrix now takes `start_adj` (the maximum likelihood solution) as the starting point to avoid inefficient sampling.

2. **Adaptive Sampling**  
   - The sampling of `zeros` and `nonzeros` is now adapted to the **newly accepted adjacency matrix** (`adj_s`) rather than the static `start_adj`.

3. **Efficiency Enhancements**  
   - Skip updates on rejection since the graph remains unchanged.  
   - Use `find_ids(adj_mat)` instead of `find_ids(start_adj)` to ensure edge modifications are correctly tracked after acceptance.

## Minor Changes:
- Example in bggm_missing.R reintroduced; Was removed due to irreproducible CRAN error.

- Addressed CRAN check error when building vignettes (removed example in bggm_missing.R -- will put it back once it's accepted to CRAN). 
- Downgraded required R version to 4.0.0

# BGGM 2.1.3 
- Replaced dprecated armadillo function `conv_to<>::from` with `as_scalar`
- `prior_sd`: Adjusted computation of delta. Also, changed default value for estimation: sqrt(1/3) resulting in delta = 2. For model testing default is more tight, at `sigma_sd` = 0.5, resulting in delta = 3. 
- `prior_sd` is now limited to range 0 -- sqrt(1/2)

# BGGM 2.1.2
- The prior_sd (or rho_sd in var_estimate() ) is limited to ranges between 0 and sqrt(1/8). These values ensure that delta does not go below 1.
- *Critical*: select() did not return partial correlations, but Fisher-z values in summary(). Fisher values are transformed back to correlation metric. This fixes #90, see [changes](https://github.com/donaldRwilliams/BGGM/commit/a264c440006069e5f171494d9618bae57f4d6566).
- Upgraded deprecated ggplot guides() argument
- Resolved non positive definite initialization matrix in wishrnd() in copula models when NA's are present in observed variables (fixes #89). See changes [here](https://github.com/donaldRwilliams/BGGM/commit/d57a5ebabd665907622a1c635ca32b5c6c913184)

# BGGM 2.1.1
BFpack dependency error fixed. 

# BGGM 2.0.1
This version of BGGM included changes based on the JOSS reviews: see [here](https://github.com/openjournals/joss-reviews/issues/2111) for 
the overview and [here](https://github.com/donaldRwilliams/BGGM/issues?q=is%3Aissue+is%3Aclosed) for specific issues.


# BGGM 2.0.0

**BGGM** was almost completely rewritten for version `2.0.0`. This was due to adding support 
for binary, ordinal, and mixed data, which required that the methods be written in `c ++`. 
Unfortunately, as a result, lots of code from version `1.0.0` is broken.

## Added features

* Full support for binary, ordinal, and mixed data. This is implemented with the argument `type`

* `roll_your_own`: compute custom network statistics from a weighted adjacency matrix or a partial 
correlation matrix

* `pcor_to_cor`: convert the sampled partial correlation matrices into correlation matrices. 

* `zero_order_cors`: compute zero order correlations 

* `convergence`: acf and trace plots

* `posterior_samples`: extract posterior samples

* `regression_summary`: summarize multivariate regression

* `pcor_sum`: Compute and compare partial correlation sums

* `weighted_adj_mat`: Extract the Weighted Adjacency Matrix

* `pcor_mat`: 	Extract the Partial Correlation Matrix

* Five additional data sets were added.

## Extensions
* `ggm_compare_ppc`: added option for custom network statistics

* Added option to control for variables with `formula`

* A progress bar was added to many functions


# BGGM 1.0.0

Initial CRAN release
