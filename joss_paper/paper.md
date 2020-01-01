---
title: 'BGGM: A R package for Gaussian Graphical Models'
tags:
  - R
  - Gaussian graphical models
  - Bayesian
  - Bayes factor
  - partial correlation
authors:
  - name: Donald R. Williams
    affiliation: "Department of Psychology, University of California, Davis" # (Multiple affiliations must be quoted)
  - name: Joris Mulder
    affiliation: "Department of Methodology and Statistics, Tilburg University"
date: 18 August 2019
bibliography: paper.bib
---

# Summary

Gaussian graphical models (GGM) allow for learning conditional dependence structures that are encoded by non-zero partial correlations. Whereas there are several R packages for classical (i.e., frequentist) methods, there are only two that implement a Bayesian approach. These are exclusively focused on identifying the graphical structure. Thepackage **BGGM** not only contains novel Bayesian methods for this purpose, but it also includes Bayesian methodology for extending inference beyond identifying non-zero relations. **BGGM** is built around two Bayesian approaches for inference--estimation and hypothesis testing. The former focuses on the posterior distribution and includes extensions to assess predictability, as well as methodology to compare partial correlations. The latter includes methods for Bayesian hypothesis testing, in both exploratory and confirmatory contexts, with the novel matrix-$F$ prior distribution. This allows for testing order and equality constrained hypotheses, as well as a combination of both with the Bayes factor. Further, there are several approaches for comparing GGMs across any number of groups. 

**BGGM** is built around two general approach for inference--estimation and hypothesis testing. The former makes use of the posterior distribution, whereas the later employs the Bayes factor. We first describe the estimation methods and then the hypothesis testing methods. We then describe diagonistic checks for the fitted models.

## Estimation
### Structure Learning
### Predictability
### Comparing GGMs

## Hypothesis Testing
### Structure Learning
### Confirmatory Hypothesis Testing
### Comparing GGMs

## Model Checking
### Regression Diagnostics
### Posterior Predictive Checks
### MCMC Convergence


The following methods are described in @Williams2018bayes. An analytic solution for determing $E$ is implemented with 

```r
# data (p = 5 for demonstrative purposes)
> Y <- BGGM::bfi[,1:5]
# fit model
> fit_analytic <- estimate(Y, analytic = T)
# select edge set
> E <- select(fit_analytic, ci_width = 0.95)
# summary
> summary(E)
# output

BGGM: Bayesian Gaussian Graphical Models 
--- 
Type: Selected Graph (Analytic Solution) 
Credible Interval: 95 % 
Connectivity: 80 % 
--- 
Call:
select.estimate(x = fit_analytic, ci_width = 0.95)
--- 
Selected:
 
Partial correlations 
 
      1     2     3    4    5
1  0.00 -0.24 -0.11 0.00 0.00
2 -0.24  0.00  0.29 0.16 0.16
3 -0.11  0.29  0.00 0.18 0.36
4  0.00  0.16  0.18 0.00 0.12
5  0.00  0.16  0.36 0.12 0.00
--- 
 
Adjacency 
 
  1 2 3 4 5
1 0 1 1 0 0
2 1 0 1 1 1
3 1 1 0 1 1
4 0 1 1 0 1
5 0 1 1 1 0
--- 
```
