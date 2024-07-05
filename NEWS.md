# BGGM 2.1.3
- `prior_sd`: Default value for estimation: sqrt(1/3) resulting in delta = 2. For model testing default is more tight, at `sigma_sd` = 0.5, resulting in delta = 3. 
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
