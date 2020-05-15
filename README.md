
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
2018) and (2) hypothesis testing (Williams and Mulder 2019). The key
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
algorithims being written in `c++` and support for all data types. The
developmental version is essentially **BGGM** version 2.0.0.

# Overview

The methods in **BGGM** build upon existing algorithms that are
well-known in the literature. The central contribution of **BGGM** is to
extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution
    (Mulder and Pericchi 2018)
    
      - [Estimation](#bayesian-estimation) (Williams 2018)

2.  Bayesian hypothesis testing with the matrix-F prior distribution
    (Mulder and Pericchi 2018).
    
      - [Exploratory hypothesis testing](#Exploratory) (Williams and
        Mulder 2019).
    
      - [Confirmatory hypothesis testing](#Confirmatory) (Williams and
        Mulder 2019).

3.  Comparing Gaussian graphical models
    
      - [Partial correlation
        differences](#partial-correlation-differences) (Williams 2018)
    
      - [Posterior predictive check](#posterior-predictive%20check)
        (Williams et al. 2020)
    
      - [Exploratory hypothesis testing](#exploratory-\(groups\))
        (Williams et al. 2020)
    
      - [Confirmatory hypothesis testing](#confirmatory-\(groups\))
        (Williams et al. 2020)

4.  Extending inference beyond the conditional (in)dependence structure
    
      - [Predictability](#Predictability) (Williams 2018)
    
      - [Posterior uncertainty intervals](#posterior-uncertainty) in the
        partial correlations (Williams 2018)
    
      - [Custom Network Statistics](#custom-network-statistics)

The computationally intensive tasks are written in `c++` via the `R`
package **Rcpp** (Eddelbuettel et al. 2011) and the `c++` library
**Armadillo** (Sanderson and Curtin 2016). The Bayes factors are
computed with the `R` package **BFpack** (Mulder et al. 2019).
Furthermore, there are plotting functions for each method, control
variables can be included in the model (e.g., `~ gender`), and there is
support for missing values (see `bggm_missing`).

## Supported Data Types

  - **Continuous**: the continuous method was described in Williams
    (2018).

  - **Binary**: the binary method builds directly upon (2012) that, in
    turn, built upon the approaches of \[-lawrence2008bayesian\] and
    \[-webb2008bayesian;textual\] (to name a few).

  - **Ordinal**: the ordinal methods require sampling thresholds. There
    are two approach included in **BGGM**. The customary approach
    described in (1993) (the default) andthe ‘Cowles’ algorithm
    described in (- Cowles 1996).

  - **Mixed**: the mixed data (a combination of discrete and continuous)
    method was introduced in \[(2007);textual\]. This is a
    semi-parametric copula model (i.e., a copula GGM) based on the
    ranked likelihood. Note that this can be used for data consisting
    entirely of ordinal data (not restricted to “mixed” data).

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

<div id="ref-albert1993bayesian">

Albert, James H, and Siddhartha Chib. 1993. “Bayesian Analysis of Binary
and Polychotomous Response Data.” *Journal of the American Statistical
Association* 88 (422): 669–79.

</div>

<div id="ref-cowles1996accelerating">

Cowles, Mary Kathryn. 1996. “Accelerating Monte Carlo Markov Chain
Convergence for Cumulative-Link Generalized Linear Models.” *Statistics
and Computing* 6 (2): 101–11.

</div>

<div id="ref-eddelbuettel2011rcpp">

Eddelbuettel, Dirk, Romain François, J Allaire, Kevin Ushey, Qiang Kou,
N Russel, John Chambers, and D Bates. 2011. “Rcpp: Seamless R and C++
Integration.” *Journal of Statistical Software* 40 (8): 1–18.

</div>

<div id="ref-hoff2007extending">

Hoff, Peter D. 2007. “Extending the Rank Likelihood for Semiparametric
Copula Estimation.” *The Annals of Applied Statistics* 1 (1): 265–83.

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

<div id="ref-talhouk2012efficient">

Talhouk, Aline, Arnaud Doucet, and Kevin Murphy. 2012. “Efficient
Bayesian Inference for Multivariate Probit Models with Sparse Inverse
Correlation Matrices.” *Journal of Computational and Graphical
Statistics* 21 (3): 739–57.

</div>

<div id="ref-Williams2019">

Williams, Donald R. 2018. “Bayesian Estimation for Gaussian Graphical
Models: Structure Learning, Predictability, and Network Comparisons.”
*arXiv*. <https://doi.org/10.31234/OSF.IO/X8DPR>.

</div>

<div id="ref-Williams2019_bf">

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
