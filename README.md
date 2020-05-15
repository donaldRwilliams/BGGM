
<img src="joss_paper/hex.png" width = 250 />

# BGGM: Bayesian Gaussian Graphical Models

[![CRAN
Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build
Status](https://travis-ci.org/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.org/donaldRwilliams/BGGM)

The `R` package **BGGM** provides tools for making Bayesian inference in
Gaussian graphical models (GGM). The methods are organized around two
general approaches for Bayesian inference: (1) estimation (Williams
2019) and (2) hypothesis testing (Williams and Mulder 2019). The key
distinction is that the former focuses on either the posterior or
posterior predictive distribution, whereas the latter focuses on model
comparison with the Bayes factor.

## Installation

To install the latest release version (1.0.0) from CRAN use

``` r
install.packages("BGGM")    
```

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("donaldRwilliams/BGGM")
```

Note that the developmental version is recommended, due to the
algorithims being written in `c++` and support for all data types (e.g.,
ordinal). As of 5/15/2020, the developmental version is essentially
**BGGM** version 2.0.0.

# Overview

The methods in **BGGM** build upon existing algorithms that are
well-known in the literature. The central contribution of **BGGM** is to
extend those approaches:

<!-- and (2) hypothesis testing -->

<!-- \insertCite{Williams2019_bf}{BGGM}. The key distinction is that the former focuses on either the -->

<!-- posterior or posterior predictive distribution, whereas the latter focuses on model comparison -->

<!-- with the Bayes factor. The methods in \strong{BGGM} build upon existing algorithms that are well-known in the literature. The central contribution of \strong{BGGM} is to extend those approaches: -->

## References

<div id="refs" class="references">

<div id="ref-Williams2019">

Williams, Donald R. 2019. “Bayesian Estimation for Gaussian Graphical
Models: Structure Learning, Predictability, and Network Comparisons.”
*arXiv*. <https://doi.org/10.31234/OSF.IO/X8DPR>.

</div>

<div id="ref-Williams2019_bf">

Williams, Donald R., and Joris Mulder. 2019. “Bayesian Hypothesis
Testing for Gaussian Graphical Models: Conditional Independence and
Order Constraints.” *PsyArXiv*. <https://doi.org/10.31234/osf.io/ypxd8>.

</div>

</div>
