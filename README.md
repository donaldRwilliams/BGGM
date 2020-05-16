
<img src="joss_paper/hex.png" width = 250 />

# BGGM: Bayesian Gaussian Graphical Models

[![CRAN
Version](http://www.r-pkg.org/badges/version/BGGM)](https://cran.r-project.org/package=BGGM)
[![Downloads](https://cranlogs.r-pkg.org/badges/BGGM)](https://cran.r-project.org/package=BGGM)
[![Build
Status](https://travis-ci.org/donaldRwilliams/BGGM.svg?branch=master)](https://travis-ci.org/donaldRwilliams/BGGM)

The `R` package **BGGM** provides tools for making Bayesian inference in
Gaussian graphical models (GGM). The methods are organized around two
general approaches for Bayesian inference: (1) estimation and (2)
hypothesis testing. The key distinction is that the former focuses on
either the posterior or posterior predictive distribution (Gelman, Meng,
and Stern 1996; see section 5 in Rubin 1984) , whereas the latter
focuses on model comparison with the Bayes factor (Jeffreys 1961; Kass
and Raftery 1995).

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
algorithms being written in `c++` and support for all data types. The
developmental version is essentially **BGGM** version 2.0.0.

## Overview

The methods in **BGGM** build upon existing algorithms that are
well-known in the literature. The central contribution of **BGGM** is to
extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution
    (Mulder and Pericchi 2018)
    
      - [Estimation](#bayesian-estimation) (Williams 2018)

2.  Bayesian hypothesis testing with the matrix-F prior distribution
    (Williams and Mulder 2019)
    
      - [Exploratory hypothesis testing](#Exploratory)
    
      - [Confirmatory hypothesis testing](#Confirmatory)

3.  Comparing Gaussian graphical models (Williams 2018; Williams et al.
    2020)
    
      - [Partial correlation
        differences](#partial-correlation-differences)
    
      - [Posterior predictive check](#posterior-predictive-check)
    
      - [Exploratory hypothesis testing](#exploratory-groups)
    
      - [Confirmatory hypothesis testing](#confirmatory-groups)

4.  Extending inference beyond the conditional (in)dependence structure
    (Williams 2018)
    
      - [Predictability](#Predictability)
    
      - [Posterior uncertainty intervals](#posterior-uncertatiny) in the
        partial correlations
    
      - [Custom Network Statistics](#custom-network-statistics)

The computationally intensive tasks are written in `c++` via the `R`
package **Rcpp** (Eddelbuettel et al. 2011) and the `c++` library
**Armadillo** (Sanderson and Curtin 2016). The Bayes factors are
computed with the `R` package **BFpack** (Mulder et al. 2019).
Furthermore, there are [plotting](#example-network-plot) functions for
each method, control variables can be included in the model (e.g., `~
gender`), and there is support for missing values (see `bggm_missing`).

## Supported Data Types

  - **Continuous**: The continuous method was described in Williams
    (2018). Note that this is based on the customary [Wishart
    distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

  - **Binary**: The binary method builds directly upon Talhouk, Doucet,
    and Murphy (2012) that, in turn, built upon the approaches of
    Lawrence et al. (2008) and Webb and Forster (2008) (to name a few).

  - **Ordinal**: The ordinal methods require sampling thresholds. There
    are two approach included in **BGGM**. The customary approach
    described in Albert and Chib (1993) (the default) and the ‘Cowles’
    algorithm described in Cowles (1996).

  - **Mixed**: The mixed data (a combination of discrete and continuous)
    method was introduced in Hoff (2007). This is a semi-parametric
    copula model (i.e., a copula GGM) based on the ranked likelihood.
    Note that this can be used for *only* ordinal data (not restricted
    to “mixed” data).

## Illustrative Examples

The following includes brief examples for *some* of the methods in
**BGGM**.

### Bayesian Estimation

#### Posterior Sampling

An ordinal GGM is estimated with

``` r
# data
Y <- ptsd[,1:5] + 1

# ordinal
fit <- estimate(Y, type = "ordinal", 
                analytic = FALSE)
```

Notice the `+ 1`. This is required, because the first category must be
`1` when `type = "ordinal"`. The partial correlations can the be
summarized with

``` r
summary(fit)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 250 
#> Observations (n):
#> Nodes (p): 5 
#> Relations: 10 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "ordinal", analytic = FALSE, iter = 250)
#> --- 
#> Estimates:
#>  Relation Post.mean Post.sd Cred.lb Cred.ub
#>    B1--B2     0.258   0.079   0.105   0.418
#>    B1--B3     0.028   0.086  -0.127   0.189
#>    B2--B3     0.517   0.058   0.406   0.616
#>    B1--B4     0.356   0.070   0.210   0.486
#>    B2--B4    -0.076   0.078  -0.219   0.063
#>    B3--B4     0.246   0.077   0.107   0.385
#>    B1--B5     0.131   0.080  -0.020   0.279
#>    B2--B5     0.127   0.083  -0.040   0.284
#>    B3--B5     0.202   0.079   0.063   0.366
#>    B4--B5     0.349   0.070   0.209   0.474
#> ---
```

The returned object can also be plotted, which allows for visualzing the
posterior uncertainty interval for each relation. This summary can the
be plotted. An example is provided below in [Posterior uncertainty
intervals](#posterior-uncertatiny). The partial correlation matrix is
accesed with

``` r
pcor_mat(fit)
```

|    |    B1 |      B2 |    B3 |      B4 |    B5 |
| -- | ----: | ------: | ----: | ------: | ----: |
| B1 | 0.000 |   0.258 | 0.028 |   0.356 | 0.131 |
| B2 | 0.258 |   0.000 | 0.517 | \-0.076 | 0.127 |
| B3 | 0.028 |   0.517 | 0.000 |   0.246 | 0.202 |
| B4 | 0.356 | \-0.076 | 0.246 |   0.000 | 0.349 |
| B5 | 0.131 |   0.127 | 0.202 |   0.349 | 0.000 |

The graph is selected with

``` r
select(fit)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Analytic: FALSE 
#> Formula:  
#> Posterior Samples: 250 
#> Credible Interval: 95 % 
#> --- 
#> Call: 
#> estimate(Y = Y, type = "ordinal", analytic = FALSE, iter = 250)
#> --- 
#> Selected:
#> 
#>       1     2     3     4     5
#> 1 0.000 0.258 0.000 0.356 0.000
#> 2 0.258 0.000 0.517 0.000 0.000
#> 3 0.000 0.517 0.000 0.246 0.202
#> 4 0.356 0.000 0.246 0.000 0.349
#> 5 0.000 0.000 0.202 0.349 0.000
#> ---
```

and then plotted

``` r
# "communities"
comm <- substring(colnames(Y), 1, 1)

plot(select(fit), 
     groups = comm,
     edge_magnify = 5, 
     palette = "Pastel1", 
     node_size = 12)
#> $plt
```

<img src="joss_paperunnamed-chunk-9-1.png" width="50%" style="display: block; margin: auto;" />

This basic “workflow” can be used with all methods and data types. A
more involved network plot is provided below.

#### Analytic

There is also an analytic solution that is based on the Wishart
distribution. This simple solution provides competivie performance with
“state-of-the-art” methods, assuming that *n* (observations) \> *p*
(variables) The one caveat is that it works only for `type =
"continuous"` (the default).

``` r
# analytic
fit <- estimate(Y, analytic = TRUE)

# network plot
plot(select(fit))
```

This is quite handy when (1) only the conditional dependence structure
is of interest and (2) an immediate solution is desirable. An example of
(2) is provided in [Posterior Predictive
Check](#posterior-predictive-check).

### Bayesian Hypothesis Testing

#### Exploratory

#### Confirmatory

### Comparing Gaussian Graphical Models

#### Partial Correlation Differences

#### Posterior Predictive Check

#### Exploratory (groups)

#### Confirmatory (groups)

### Beyond the Conditional (in)dependence Structure

#### Predictability

#### Posterior Uncertainty

#### Custom Network Statistics

### Example Network Plot

## Additional Features

The primary focus of **BGGM** is Gaussian graphical modeling (the
inverse covariance matrix). The residue is a suite of useful methods not
explicitly for GGMs. For example,

### Bivariate Correlations

Bivariate correlations for `binary` (tetrachoric), `ordinal`
(polychoric), `mixed` (rank based), and `continuous` (Pearson’s) data.

Here is an example for computing tetrachoric correlations:

``` r
# binary data
Y <- women_math[1:500,]

cors <- zero_order_cors(Y, type = "binary")

cors$R
```

|   |       1 |       2 |       3 |       4 |       5 |       6 |
| - | ------: | ------: | ------: | ------: | ------: | ------: |
| 1 |   1.000 | \-0.198 |   0.506 |   0.122 | \-0.140 |   0.098 |
| 2 | \-0.198 |   1.000 | \-0.482 | \-0.013 | \-0.146 | \-0.146 |
| 3 |   0.506 | \-0.482 |   1.000 |   0.310 | \-0.343 |   0.351 |
| 4 |   0.122 | \-0.013 |   0.310 |   1.000 | \-0.363 |   0.169 |
| 5 | \-0.140 | \-0.146 | \-0.343 | \-0.363 |   1.000 | \-0.194 |
| 6 |   0.098 | \-0.146 |   0.351 |   0.169 | \-0.194 |   1.000 |

The object `cors` also includes the sampled correlation matrices (in
this case 250) in an array.

### Multivariate Regression

Multivariate regression for binary (probit), ordinal (probit), mixed
(rank likelihood), and continuous data.

Here is an example for a multivariate probit model with an ordinal
outcome, where `E5` (“take charge”) and `N5` (“panic easily”) are
predicted by `gender` and `education`:

``` r
# personality data
Y <- bfi

# variables
Y <- subset(Y, select = c("E5", "N5", 
                          "gender", "education"))


mv_probit <- estimate(Y, formula = ~ gender + as.factor(education), 
                      type = "ordinal")
```

Note that **BGGM** does not use the customary `model.matrix`
formulation. This is for good reason, as each variable in the GGM does
not need to be written out. Here we effectively “tricked” **BGGM** to
fit a multivariate probit model (each variable included in `formula` is
removed from `Y`).

``` r
regression_summary(mv_probit)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: ordinal 
#> Formula: ~ gender + as.factor(education) 
#> --- 
#> Coefficients: 
#>  
#> E5 
#>                       Post.mean Post.sd Cred.lb Cred.ub
#> (Intercept)               1.852   0.533   1.049   3.142
#> gender                    0.169   0.066   0.065   0.295
#> as.factor(education)2     0.215   0.109   0.024   0.437
#> as.factor(education)3     0.271   0.104   0.089   0.445
#> as.factor(education)4     0.206   0.103   0.019   0.404
#> as.factor(education)5     0.345   0.128   0.120   0.593
#> --- 
#> N5 
#>                       Post.mean Post.sd Cred.lb Cred.ub
#> (Intercept)               0.210   0.114  -0.012   0.434
#> gender                    0.502   0.140   0.291   0.835
#> as.factor(education)2    -0.127   0.103  -0.345   0.058
#> as.factor(education)3    -0.104   0.081  -0.258   0.034
#> as.factor(education)4    -0.218   0.104  -0.427  -0.024
#> as.factor(education)5    -0.229   0.103  -0.449  -0.038
#> --- 
#> Residual Correlation Matrix: 
#>       E5    N5
#> E5  1.00 -0.18
#> N5 -0.18  1.00
#> ---
```

This basic idea can also be used to fit regression models with a single
outcome.

## Note on Conditional (In)dependence Models for Latent Data

All of the data types (besides continuous) model latent data. That is,
unobserved data that is assumed to be Gaussian distributed. For example,
a tetrachoric correlation (binary data) is a special case of a
polychoric correlation (ordinal data). Both relations are between
"theorized normally distributed continuous *latent* variables
[Wikepedia](https://en.wikipedia.org/wiki/Polychoric_correlation). In
both instances, the corresponding partial correlation between observed
variables is conditioned on the remaining variables in the *latent*
space. This implies that interpretation is similar to continuous data,
but with respect to latent variables. We refer interested users to (see
page 2364, section 2.2, in Webb and Forster 2008).

## High Dimensional Data?

**BGGM** was built specifically for social-behavioral scientists. Of
course, the methods can be used by all researchers. However, there is
currently *not* support for high-dimensional data (i.e., more variables
than observations) that are common place in, say, the genetics
literature. These data are rare in the social-behavioral sciences. In
the future, support for high-dimensional data may be added to **BGGM**.

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

<div id="ref-Gelman1996a">

Gelman, Andrew, Xiao-Li Meng, and Hal Stern. 1996. “Posterior predictive
assessment of model fitness via realized discrepancies. Vol.6, No.4.”
*Statistica Sinica* 6 (4): 733–807. <https://doi.org/10.1.1.142.9951>.

</div>

<div id="ref-hoff2007extending">

Hoff, Peter D. 2007. “Extending the Rank Likelihood for Semiparametric
Copula Estimation.” *The Annals of Applied Statistics* 1 (1): 265–83.

</div>

<div id="ref-Jeffreys1961">

Jeffreys, Harold. 1961. *The theory of probability*. Oxford: Oxford
University Press.

</div>

<div id="ref-Kass1995">

Kass, Robert E, and Adrian E Raftery. 1995. “Bayes Factors.” *Journal of
the American Statistical Association* 90 (430): 773–95.

</div>

<div id="ref-lawrence2008bayesian">

Lawrence, Earl, Derek Bingham, Chuanhai Liu, and Vijayan N Nair. 2008.
“Bayesian Inference for Multivariate Ordinal Data Using Parameter
Expansion.” *Technometrics* 50 (2): 182–91.

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

<div id="ref-rubin1984bayesianly">

Rubin, Donald B. 1984. “Bayesianly Justifiable and Relevant Frequency
Calculations for the Applied Statistician.” *The Annals of Statistics*,
1151–72.

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

<div id="ref-webb2008bayesian">

Webb, Emily L, and Jonathan J Forster. 2008. “Bayesian Model
Determination for Multivariate Ordinal and Binary Data.” *Computational
Statistics & Data Analysis* 52 (5): 2632–49.

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
