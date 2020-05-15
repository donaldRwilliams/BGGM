
<img src="joss_paper/hex.png" width = 250 />

# BGGM: Bayesian Gaussian Graphical Models

[![CRAN
Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build
Status](https://travis-ci.org/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.org/donaldRwilliams/BGGM)

The `R` package **BGGM** provides tools for making Bayesian inference in
Gaussian graphical models (GGM). The methods are organized around two
general approaches for Bayesian inference: (1) estimation (D. R.
Williams 2019) and (2) hypothesis testing
(<span class="citeproc-not-found" data-reference-id="Williams201b">**???**</span>).
The key distinction is that the former focuses on either the posterior
or posterior predictive distribution, whereas the latter focuses on
model comparison with the Bayes factor.

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

1.  Bayesian estimation with the novel matrix-F prior distribution
    (Mulder and Pericchi 2018)
    
      - [Estimation](#bayesian-estimation) (D. R. Williams 2019)

2.  Bayesian hypothesis testing with the matrix-F prior distribution
    (Mulder and Pericchi 2018).
    
      - [Exploratory hypothesis testing](#Exploratory) (D. R. Williams
        and Mulder 2019).
    
      - [Confirmatory hypothesis testing](#Confirmatory) (D. R. Williams
        and Mulder 2019).

3.  Comparing Gaussian graphical models
    
      - [Partial correlation
        differences](#partial-correlation-differences) (D. R. Williams
        2019)
    
      - [Posterior predictive check](#posterior-predictive%20check)
        (Williams et al. 2020)
    
      - [Exploratory hypothesis testing](#exploratory-\(groups\))
        (Williams et al. 2020)
    
      - [Confirmatory hypothesis testing](#confirmatory-\(groups\))
        (Williams et al. 2020)

4.  Extending inference beyond the conditional (in)dependence structure
    
      - [Predictability](#Predictability) (D. R. Williams 2019)
    
      - [Posterior uncertainty intervals](#posterior-uncertainty) in the
        partial correlations (D. R. Williams 2019)
    
      - [Custom Network Statistics](#custom-network-statistics)

The computationally intensive tasks are written in `c++` via the `R`
package **Rcpp** (Eddelbuettel et al. 2011) and the `c++` library
**Armadillo** (Sanderson and Curtin 2016). The Bayes factors are
computed with the `R` package **BFpack** (Mulder et al. 2019).
Furthermore, there are plotting functions for each method, control
variables can be included in the model (e.g., `~ gender`), and there is
support for missing values (see `bggm_missing`).

# Illustrative Examples

## Bayesian Estimation

### Posterior Sampling

### Analytic

## Bayesian Hypothesis Testing

### Exploratory

### Confirmatory

## Comparing Gaussian Graphical Models

### Partial Correlation Differences

### Posterior Predictive Check

### Exploratory (groups)

### Confirmatory (groups)

## Beyond the Conditional (in)dependence Structure

### Predictability

### Posterior Uncertatiny

### Custom Network Statistics

## References

<div id="refs" class="references">

<div id="ref-eddelbuettel2011rcpp">

Eddelbuettel, Dirk, Romain François, J Allaire, Kevin Ushey, Qiang Kou,
N Russel, John Chambers, and D Bates. 2011. “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software* 40 (8): 1–18.

</div>

<div id="ref-mulder2019bfpack">

Mulder, Joris, Xin Gu, Anton Olsson-Collentine, Andrew Tomarken, Florian
Böing-Messing, Herbert Hoijtink, Marlyne Meijerink, et al. 2019.
“BFpack: Flexible Bayes Factor Testing of Scientific Theories in R.”
*arXiv Preprint arXiv:1911.07728*.

</div>

<div id="ref-Mulder2018">

Mulder, Joris, and Luis Pericchi. 2018. “The Matrix-F Prior for
Estimating and Testing Covariance Matrices.” *Bayesian Analysis*, no. 4:
1–22. <https://doi.org/10.1214/17-BA1092>.

</div>

<div id="ref-sanderson2016armadillo">

Sanderson, Conrad, and Ryan Curtin. 2016. “Armadillo: A Template-Based
C++ Library for Linear Algebra.” *Journal of Open Source Software* 1
(2): 26.

</div>

<div id="ref-Williams2019a">

Williams, Donald R. 2019. “Bayesian Estimation for Gaussian Graphical
Models: Structure Learning, Predictability, and Network Comparisons.”
*arXiv*. <https://doi.org/10.31234/OSF.IO/X8DPR>.

</div>

<div id="ref-Williams2019b">

Williams, Donald R., and Joris Mulder. 2019. “Bayesian Hypothesis
Testing for Gaussian Graphical Models: Conditional Independence and
Order Constraints.” *PsyArXiv*. <https://doi.org/10.31234/osf.io/ypxd8>.

</div>

<div id="ref-williams2020comparing">

Williams, Donald R., Philippe Rast, Luis R Pericchi, and Joris Mulder.
2020. “Comparing Gaussian Graphical Models with the Posterior Predictive
Distribution and Bayesian Model Selection.” *Psychological Methods*.

</div>

</div>
