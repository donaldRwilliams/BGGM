---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'BGGM: Bayesian Gaussian Graphical Models in R'
tags:
- Gaussian graphical models
- Bayesian
- Bayes factor
- partial correlation
- R
authors:
  - name: Donald R. Williams
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Joris Mulder
    affiliation: 2
affiliations:
 - name: Department of Psychology, University of California, Davis
   index: 1
 - name: Department of Methodology and Statistics, Tilburg University
   index: 2
citation_author: Williams and Mulder
date: 05 May 2020
year: 2020
bibliography: inst/REFERENCES.bib
---

# BGGM: Bayesian Gaussian Graphical Models
The `R` package **BGGM** provides tools for making Bayesian inference in 
Gaussian graphical models (GGM). The methods are organized around two general 
approaches for Bayesian inference: (1) estimation and (2) hypothesis 
testing. The key distinction is that the former focuses on either 
the posterior or posterior predictive distribution [@Gelman1996a; see section 
5 in @rubin1984bayesianly] , whereas the latter focuses on model comparison with the Bayes factor [@Jeffreys1961; @Kass1995].

# Overview
The methods in **BGGM** build upon existing algorithms that are well-known in the literature.
The central contribution of **BGGM** is to extend those approaches:

1.  Bayesian estimation with the novel matrix-F prior distribution [@Mulder2018]
  
    + [Estimation](https://github.com/donaldRwilliams/BGGM#bayesian-estimation) [@Williams2019]

2. Bayesian hypothesis testing with the matrix-F prior distribution [@Williams2019_bf]

    + [Exploratory hypothesis testing](https://github.com/donaldRwilliams/BGGM#Exploratory)
  
    + [Confirmatory hypothesis testing](https://github.com/donaldRwilliams/BGGM#Confirmatory)
    
3. Comparing Gaussian graphical models [@Williams2019; @williams2020comparing]
    
    + [Partial correlation differences](https://github.com/donaldRwilliams/BGGM#partial-correlation-differences) 
    
    + [Posterior predictive check](https://github.com/donaldRwilliams/BGGM#posterior-predictive-check)
    
    + [Exploratory hypothesis testing](https://github.com/donaldRwilliams/BGGM#exploratory-groups) 
    
    + [Confirmatory hypothesis testing](https://github.com/donaldRwilliams/BGGM#confirmatory-groups)

4. Extending inference beyond the conditional (in)dependence structure [@Williams2019]

    +  [Predictability](https://github.com/donaldRwilliams/BGGM#Predictability) 
    
    +  [Posterior uncertaintyintervals](https://github.com/donaldRwilliams/BGGM#partial-correlation-differences) for the 
       partial correlations
       
    +  [Custom Network Statistics](https://github.com/donaldRwilliams/BGGM#custom-network-statistics)
    
    
## Supported Data Types

* **Continuous**: The continuous method was described in  @Williams2019. Note that 
                  this is based on the customary[Wishartdistribution](https://en.wikipedia.org/wiki/Wishart_distribution).

* **Binary**: The binary method builds directly upon @talhouk2012efficient
  that, in turn, built upon the approaches of @lawrence2008bayesian and
  @webb2008bayesian (to name a few).
  
* **Ordinal**: The ordinal methods require sampling thresholds. There are two approach 
   included in **BGGM**. The customary approach described in @albert1993bayesian 
   (the default) and the 'Cowles' algorithm described in @cowles1996accelerating.
   
* **Mixed**: The mixed data (a combination of discrete and continuous) method was introduced
 in @hoff2007extending. This is a semi-parametric copula model
 (i.e., a copula GGM) based on the ranked likelihood. Note that this can be used for 
 *only* ordinal data (not restricted to "mixed" data).

The computationally intensive tasks are written in `c++` via the `R` package **Rcpp** [@eddelbuettel2011rcpp] and the `c++` library **Armadillo** [@sanderson2016armadillo]. The Bayes factors are computed with the `R` package **BFpack** [@mulder2019bfpack]. Furthermore, there are [plotting](https://github.com/donaldRwilliams/BGGM#example-network-plot) functions
for each method, control variables can be included in the model (e.g., `~ gender`), 
and there is support for missing values (see `bggm_missing`).

# Acknowledgements
DRW was supported by a National Science Foundation Graduate Research Fellowship
under Grant No. 1650042 and JM was supported by a ERC Starting Grant (758791).

# References
