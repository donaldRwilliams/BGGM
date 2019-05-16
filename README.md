
<!-- README.md is generated from README.Rmd. Please edit that file -->
BGGM: Bayesian Gaussian Graphical Models
========================================

This package includes methods introduced in Williams, Rast, Pericchi, and Mulder (2019), Williams and Mulder (2019), and Williams (2018). The package is built around two Bayesian approaches for inference: estimation and hypothesis testing.

The estimation based methods are described in Williams (2018). They offer advantages compared to classical methods, in that a measure of uncertainty is provided for all parameters. For example, each node has a distribution for variance explained (i.e., Bayesian *R*<sup>2</sup>). Measures of out-of-sample prediction error are available, which also have a measure of uncertainty. The model is selected with credible interval exclusion of zero or a region of practical equivalence. This allows for computing the posterior probability for an assumed *null* region--i.e., conditional independence. It is also possible to compare partial correlations.

The hypothesis testing based methods are described in Williams and Mulder (2019), and allow for testing edges (i.e., partial correlations) with the Bayes factor. One-sided hypothesis testing is also possible. These methods provide (relative) evidence for the null hypothesis. There are extensions for **confirmatory hypothesis testing** in GGMs--e.g., inequality or equality constraints on the partial correlations. This allows for comparing theoretically informed models with Bayesian model selection.

Further, it is possible to assess differences as well as similarities (i.e., the null hypothesis) between GGMs. These method were introduced in Williams, Rast, Pericchi, and Mulder (2019). Graphs are compared either with the posterior predictive distribution or Bayesian model selection. The latter allows for testing hypothesized changes in graphical structures between, for example, control and treatment groups. The posterior predictive approach is based on KL-divergence. It allows for testing the assumed (null) model of group equality for the entire graph or specific variables. These methods can be used to compare any number of GGMs.

Williams, D. R. (2018). Bayesian Inference for Gaussian Graphical Models: Structure Learning, Explanation, and Prediction. ([pre-print](https://doi.org/10.31234/osf.io/x8dpr))

Williams, D. R., & Mulder, J. (2019). Bayesian Hypothesis Testing for Gaussian Graphical Models:Conditional Independence and Order Constraints. ([pre-print](https://doi.org/10.31234/osf.io/ypxd8))

Williams, D. R., Rast, P., Pericchi, L. R., & Mulder, J. (2019). Comparing Gaussian Graphical Models with the Posterior Predictive Distribution and Bayesian Model Selection. ([pre-print](https://psyarxiv.com/yt386/))

Outline
-------

This README is organized as follows:

-   [Installation](https://github.com/donaldRwilliams/BGGM#installation)

-   [Estimation](https://github.com/donaldRwilliams/BGGM#estimation)

    -   [Structure Learning (Estimation)](https://github.com/donaldRwilliams/BGGM#structure-learning-estimation)

    -   [Edge (partial correlation) Differences](https://github.com/donaldRwilliams/BGGM#edge-differences)

    -   [Prediction](https://github.com/donaldRwilliams/BGGM#prediction)

        -   [Bayesian *R*<sup>2</sup>](https://github.com/donaldRwilliams/BGGM#bayesian-r2)

        -   [Leave-One-Out Cross-Validation](https://github.com/donaldRwilliams/BGGM#leave-one-out-cross-validation)

-   [Hypothesis Testing](https://github.com/donaldRwilliams/BGGM#hypothesis-testing)

    -   [Structure Learning (Bayes Factor)](https://github.com/donaldRwilliams/BGGM#structure-learning-bayes-factor)

        -   [Visualizing Scientific Expectations](https://github.com/donaldRwilliams/BGGM#visualizing-scientific-expectations)

        -   [Two-Sided Hypothesis Testing](https://github.com/donaldRwilliams/BGGM#two-sided-testing)

        -   [One-Sided Hypothesis Testing](https://github.com/donaldRwilliams/BGGM#one-sided-testing)

        -   [Exhaustive Hypothesis Testing](https://github.com/donaldRwilliams/BGGM#exhaustive-hypothesis-testing)

    -   Confirmatory Hypothesis Testing

        -   Order Constraints

        -   Equality Constraints

-   Comparing GGMs

    -   Posterior Predictive KL-divergence

    -   Bayesian Model Selection

Installation
============

You can install BGGM from git hub with:

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/BGGM")
```

Estimation
==========

The following methods are described in Williams (2018). They build upon two basic ideas. First, the Wishart distribution is a conjugate prior distribution for the precision matrix (inverse of the covariance matrix). This provides an analytic solution for selecting the graph, and allows for conveniently drawing posterior samples. Second, there is an exact relationship between estimating the precision matrix directly and with multiple regression. Here the individual elements, from joint posterior distribution of the precision matrix, can be converted to their respective regression counterparts. This allows for assessing nodewise (for each variable in the model) predictability.

Structure Learning (Estimation)
-------------------------------

By structure learning we are referring to selecting the graph (i.e., the edge set *E*), which consists of those edges determined to be non-zero. For demonstrative purposes, we consider a relatively small number of variables (*p* = 5).

The package **BGGM** offers a convenient analytic solution for estimating GGMs. It is implemented with:

``` r
library(BGGM)
library(ggplot2)
library(ggraph)
library(foreach)


# p = 5
Y <- BGGM::bfi[,1:5]

# analytic solution
fit_analytic <- estimate(Y, analytic = T)

# summary
summary(fit_analytic)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Estimation (Analytic Solution) 
#> Posterior Samples: 
#> Observations (n): 2709 
#> Variables (p): 5 
#> Edges: 10 
#> --- 
#> Call: 
#> estimate.default(x = Y, analytic = T)
#> --- 
#> Date: Thu May 16 15:51:01 2019
```

Note `summary(.)` provides information about the fitted model, including that the analytic solution was used, the number of observations (*n*) and variables (*p*), and the number of edges.

The edge set is then selected with:

``` r
# select the graph (edge set E)
E <- select(fit_analytic, ci_width = 0.95)

# summary of E
summary(E)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Selected Graph (Analytic Solution) 
#> Credible Interval: 95 % 
#> Connectivity: 80 % 
#> --- 
#> Call:
#> select.estimate(x = fit_analytic, ci_width = 0.95)
#> --- 
#> Selected:
#>  
#> Partial correlations 
#>  
#>       1     2     3    4    5
#> 1  0.00 -0.24 -0.11 0.00 0.00
#> 2 -0.24  0.00  0.29 0.16 0.16
#> 3 -0.11  0.29  0.00 0.18 0.36
#> 4  0.00  0.16  0.18 0.00 0.12
#> 5  0.00  0.16  0.36 0.12 0.00
#> --- 
#>  
#> Adjacency 
#>  
#>   1 2 3 4 5
#> 1 0 1 1 0 0
#> 2 1 0 1 1 1
#> 3 1 1 0 1 1
#> 4 0 1 1 0 1
#> 5 0 1 1 1 0
#> ---
```

The analytic solution works directly with the precision matrix, and thus, there is not an option to summarize the posterior distributions. This is because the non-standardized elements are in the opposite direction (±) of the partial correlations, which in our experience, can lead to confusion. To summarize the partial correlations change `analytic = T` to `analytic = F`:

``` r
# sample from posterior
fit_sampling <- estimate(Y, analytic = F)

# select the graph
E <- select(fit_sampling, ci_width = 0.95)

# summarize partial correlations
summary(E, summarize = T, digits = 2)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Selected Graph (Sampling) 
#> Credible Interval: 95 % 
#> Connectivity: 80 % 
#> --- 
#> Call:
#> select.estimate(x = fit_sampling, ci_width = 0.95)
#> --- 
#> Estimates: 
#>  
#>  egde post_mean post_sd   2.5%  97.5%
#>  1--2   -0.2406   0.018 -0.276 -0.205
#>  1--3   -0.1071   0.019 -0.145 -0.069
#>  2--3    0.2861   0.018  0.251  0.321
#>  1--4   -0.0075   0.019 -0.044  0.029
#>  2--4    0.1642   0.019  0.128  0.200
#>  3--4    0.1778   0.018  0.141  0.214
#>  1--5   -0.0089   0.019 -0.046  0.030
#>  2--5    0.1563   0.019  0.119  0.193
#>  3--5    0.3585   0.017  0.326  0.391
#>  4--5    0.1217   0.018  0.085  0.158
#> ---
```

Note that `edge` corresponds to that particular entry in the partial correlation matrix--i.e., `1--2` is the relation between the first and second variables, respectively.

**BGGM** provides several options for plotting, with each implemented as a S3 generic. For example, the partial correlations can be plotted with:

``` r
# p = 10
Y <- BGGM::bfi[,1:10]

# sampling required
fit_sampling <- estimate(Y, analytic = F)

# plot
plot_1A <- plot(fit_sampling, 
                ci_width = 0.95, 
                width = 0.1,  
                size = 2) +
            coord_cartesian() +
            theme(axis.text.x = element_text(angle = 90))
  
plot_1A
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="60%" style="display: block; margin: auto;" />

This example nicely demonstrates how the `plot` objects can be further customized with **ggplot2**. There are two options for visualizing the selected graph. The heatmap plot is generated with:

``` r
# select the graph
E <- select(fit_sampling, ci_width = 0.95)

# heatmap plot
plot_1B <- plot(E, 
                type = "heatmap", 
                lower_tri = TRUE) +
           ggtitle("Heatmap Plot") + 
           theme(plot.title = element_text(size = 15))
plot_1B
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="60%" style="display: block; margin: auto;" />

Here `lower_tri = TRUE` controls which partial correlations are plotted. In this case, only the lower triangular elements are included in the plot. This can be changed with `lower_tri = FALSE`.

On the other hand, a “network” plot can be obtained with:

``` r
# network plot
plot_1C <- plot(E, type = "network",
                layout ='circle',
                node_outer = 8,
                node_inner = 7,
                node_text_size = 4) +
           ggtitle("Network Plot") +
           theme(plot.title = element_text(size = 15))
plot_1C
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="60%" style="display: block; margin: auto;" />

A key feature of **BGGM** is extending inference beyond identifying non-zero partial correlations. The region of practical equivalence can be used for this purpose, as it allows for determining which relations are practically zero. In this case, we follow Cohen’s guidelines, wherein 0.1 is considered a small effect.This is implemented with:

``` r
# p = 10
Y <- BGGM::bfi[,1:10]

# sample from posterior
fit_sample <- estimate(Y, samples = 5000, analytic = F)

# select the graph
E <- select(fit_sample, rope = 0.1, prob = 0.95)
#> ci_width is ignored

# summary for first 10 rows
head(E, nrow = 10, summarize = T, digits = 2)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Selected Graph (Sampling) 
#> Probability: 0.95 
#> Region of Practical Equivalence:[-0.1, 0.1]
#> Connectivity: 31.1 % 
#> --- 
#> Call:
#> select.estimate(x = fit_sample, rope = 0.1, prob = 0.95)
#> --- 
#> pr_out: post prob outside of rope 
#> pr_in: post prob inside of rope 
#> --- 
#> Estimates: 
#>  
#>  egde post_mean post_sd pr_out  pr_in
#>  1--2    -0.244   0.018   1.00 0.0000
#>  1--3    -0.106   0.020   0.61 0.3880
#>  2--3     0.287   0.018   1.00 0.0000
#>  1--4    -0.014   0.020   0.00 1.0000
#>  2--4     0.161   0.019   1.00 0.0006
#>  3--4     0.161   0.019   1.00 0.0004
#>  1--5    -0.016   0.020   0.00 1.0000
#>  2--5     0.144   0.019   0.99 0.0088
#>  3--5     0.354   0.017   1.00 0.0000
#>  4--5     0.114   0.019   0.77 0.2338
#> ---
```

The argument `prob = 0.95` requires that 95 % of the posterior density be in or out of the rope to be considered practically equivalent or different from zero. With this decision rule, as seen with `head(.)`, edges `1--4` and `1--5` are practically equivalent to zero. This inference is made possible with **BGGM**.

In this case, `plot(.)` returns two objects: (1) the selected edges; (2) those for which there is support for the null values. This is implemented with:

``` r
# network plot
plts <- plot(E, type = "network",
             layout ='circle',
             node_outer = 10,
             node_inner = 9,
             node_text_size = 6)

# practically non-zero
plot_1D <- plts$plot_nonzero +
             ggtitle("Practically Non-zero") +
             theme(plot.title = element_text(size = 15))

plot_1D
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="60%" style="display: block; margin: auto;" />

``` r
# practically zero
plot_1E <- plts$plot_zero +
              ggtitle("Practically Zero") +
              theme(plot.title = element_text(size = 15))

plot_1E
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="60%" style="display: block; margin: auto;" />

We emphasize that GGMs are often thought to capture conditionally *independent* relations--i.e., evidence for the null hypothesis of no effect, conditional on the other variables in the model. However, the dominant approach assesses conditional *dependence* (*ρ*<sub>*i**j*</sub> ≠ 0), and then sets relations to zero otherwise. **BGGM** can explicitly answer the question of conditional independence.

Edge Differences
----------------

Differences between partial correlations are often tested in GGMs; for example, with a classical (i.e., frequentist) approach that is implemented in **bootnet**. One contribution of **BGGM** is providing Bayesian analogs for commonly used methods, as well as extensions to those methods. In this case, we can use posterior probabilities to determine which edges are practically equivalent. This is implemented with:

``` r
# edge differences
edge_difference <- edge_compare(fit_sample, contrast = "all", ci_width = 0.95, rope = 0.1)

# summary for first 5 contrasts
head(edge_difference, nrow = 5)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Edge comparison(s) 
#> Credible Interval: 95 % 
#> Region of Practical Equivalence:[-0.1, 0.1]
#> --- 
#> Call:
#> edge_compare.estimate(x = fit_sample, contrast = "all", ci_width = 0.95, 
#>     rope = 0.1)
#> --- 
#> Estimates: 
#>  
#>   contrast post_mean post_sd pr_out pr_in
#>  1--2-1--3    -0.139   0.031  0.891 0.109
#>  1--2-2--3    -0.531   0.024  1.000 0.000
#>  1--2-1--4    -0.230   0.029  1.000 0.000
#>  1--2-2--4    -0.405   0.026  1.000 0.000
#>  1--2-3--4    -0.405   0.027  1.000 0.000
#> ---
```

This output includes the posterior mean and standard deviation for each difference. Further, `pr_in` is the proportion of samples between (±) 0.1. This can be interpreted as the posterior probability of practical equivalence, which has been defined with the argument `rope = 0.1`. Further, this powerful function can be used to assess specific contrasts. This can be accomplished, for example, with 5--1 - 6--10. Note that care must be taken when specifying the contrasts, as an error will arise if they are not in the proper format.

The object `edge_difference` can the be plotted with:

``` r
# plot contrasts
plot_diff <- plot(edge_difference, prob = .99)

# practically different
plot_2A <- plot_diff$plt_nonzero +
           ggtitle("Practically Different") +
           theme(axis.text.y = element_blank())
plot_2A
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="60%" style="display: block; margin: auto;" />

``` r
# practically equivalent
plot_2B <- plot_diff$plt_zero +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  ggtitle("Practically Equivalent") +
  theme(axis.text.y = element_blank())

plot_2B
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="60%" style="display: block; margin: auto;" />

This shows the central idea behind the region of practical equivalence, which is highlighted in grey. Ideally only a few contrasts would be examined in light of a guiding theory. To this end, the option `contrast = "all"` may be removed altogether from **BGGM**.

Prediction
----------

The following is based on the correspondence between the elements of the precision matrix and multiple regression. In the context of GGMs, using regression to select edges is referred to as “neighborhood” selection. On the other hand, the method described in Williams (2018) works directly with either the posterior distribution for the precision matrix or the maximum a posteriori estimates. These are then converted to the corresponding regression coefficients and residual variances. It follows that **BGGM** can also be used for the purpose of multiple regression–i.e.,

``` r
# p = 10
Y <- BGGM::bfi[,1:10]

# sample posterior
fit <- estimate(Y, samples = 5000)

# precision to regression
coefficients(fit, node = 1, ci_width = 0.95)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Inverse to Regression 
#> --- 
#> Call: 
#> BGGM:::beta_summary(x = fit, node = node, ci_width = ci_width, 
#>     samples = samples)
#> --- 
#> Estimates: 
#>  
#>  node post_mean post_sd   2.5%  97.5%
#>     2    -0.279   0.021 -0.321 -0.238
#>     3    -0.124   0.022 -0.168 -0.080
#>     4    -0.016   0.021 -0.057  0.025
#>     5    -0.017   0.021 -0.060  0.024
#>     6     0.057   0.020  0.014  0.097
#>     7     0.081   0.021  0.042  0.124
#>     8     0.044   0.020  0.005  0.082
#>     9     0.142   0.022  0.098  0.183
#>    10    -0.030   0.022 -0.071  0.012
#> ---
```

Here `node = 1` indicates which node is summarized. This correspondence allows for computing measures of prediction error (or accuracy), including Bayesian *R*<sup>2</sup> and Bayesian leave-one-out cross-validation, each of which has a measure of uncertainty. Furthermore, when a computationally convenient option is desirable, **BGGM** includes an analytic expression for prediction error. This is also known as the predicted residual sums of squares (PRESS).

### Bayesian *R*<sup>2</sup>

In-sample Bayesian *R*<sup>2</sup> is implemented with:

``` r
# training data
Y_train <- BGGM::bfi[1:100,1:10]

# fit to training data
fit_train <- estimate(Y_train, samples = 5000)

# compute Bayes R2
train_R2 <- predict(fit_train,
                    ci_width = 0.90,
                    samples = 1000,
                    measure = "R2")

# summary for first 2 rows
head(train_R2, nrow = 2)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: In-sample predictive accuracy 
#> Measure: Variance Explained (R2) 
#> --- 
#> Call:
#> predict.estimate(fit = fit_train, ci_width = 0.9, samples = 1000, 
#>     measure = "R2")
#> --- 
#> Estimates: 
#> 
#>  node post_mean    post_sd      2.5%     97.5%
#>     1 0.1648317 0.06682872 0.0444397 0.2945292
#>     2 0.2880843 0.06937293 0.1406179 0.4114294
#> ---
```

Here `ci_width = 0.90` indicates the decision rule for setting coefficients to zero, and by default, 95 % intervals are used in the summary output. Similarly, out-of-sample Bayesian *R*<sup>2</sup> is computed with:

``` r
# test data
Y_test <-  BGGM::bfi[101:2000,1:10]

# predict test data
test_R2 <- predict(fit_train, ci_width = 0.90,
                   test_data = Y_test,
                   samples = 1000, measure = "R2")
```

The work flow is completed by visualizing Bayesian *R*<sup>2</sup> for each node–i.e.,

``` r
# prior training and test error in the same plot
plt_3A <- plot(x1 = train_R2, x2 =  test_R2, order = "test")

plt_3A
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="60%" style="display: block; margin: auto;" />

Here the nodes have been ordered by which has the best out-of-sample performance. It is also possible to have each in a separate plot by leaving `x2` empty. The `predict` object can be used to assess differences in predictive accuracy with compare(.). **BGGM** also includes mean squared error (`measure = "mse"`).

### Leave-one-out cross-validation

Bayesian leave-one-out cross-validation is implemented with:

``` r
# p = 10
Y <- BGGM::bfi[1:1000,1:10]

# sample posterior
fit_sample <- estimate(Y, samples = 5000)

# Bayesian LOO
bayes_loo <- loocv(fit_sample)

# nodewise loo summary
summary(bayes_loo)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Leave-One-Out Prediction Error (Bayesian) 
#> --- 
#> Call:
#> loocv.default(x = fit_sample)
#> --- 
#> Estimates: 
#> 
#>   node      loo   loo_se
#>     1 2573.365 48.79506
#>     2 2330.739 63.43818
#>     3 2302.289 63.62531
#>     4 2466.652 51.87561
#>     5 2416.475 55.55979
#>     6 2433.823 58.47885
#>     7 2300.522 50.78019
#>     8 2391.761 51.11185
#>     9 2296.453 51.43718
#>    10 2365.595 39.89154
#> ---
```

The results are plotted with:

``` r
# plot CV error
plt_3B <- plot(bayes_loo, size = 8) +
          theme_classic() +
          ylab("Bayesian Leave-One-Out")

plt_3B
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="60%" style="display: block; margin: auto;" />

Similarly, by setting `analytic = T`, leave-one-out prediction error can be computed analytically. This is implemented with:

``` r
# p = 10
Y <- BGGM::bfi[1:1000,1:10]

# analytic solution
fit_analytic <- estimate(Y, analytic = T)

# analytic LOO (PRESS; based on point estimates)
press_loo <- loocv(fit_analytic)

# plot CV error
plt_3C <- plot(press_loo, size = 8) +
          theme_classic() +
          ylab("PRESS: Leave-One-Out") +
          scale_y_continuous(expand = c(0, 0),
          limit = c(0, 1000))

plt_3C
```

<img src="man/figures/README-unnamed-chunk-20-1.png" width="60%" style="display: block; margin: auto;" />

This highlights the difference between the leave-one-out methods, in that the Bayesian version has a measure of uncertainty (although the order is the same). For both measures of predictive *error*, a lower value indicates a more predictable node (variable).

Hypothesis Testing
==================

The following methods were introduced in Williams and Mulder (2019). That work not only presented an exploratory approach using the Bayes factor, but it also proposed methodology for confirmatory hypothesis testing in GGMs. The latter provides an alternative to data driven model selection that is commonplace in the GGM literature, and in particular, it allows for comparing theoretical models. Further, the novel matrix−*F* prior distribution is used for the partial correlations, which offers more flexibility than the Wishart distribution. This approach builds upon (Mulder2016), where the focus was on correlations (as opposed to *partial* correlations). In particular, **BGGM** allows for Bayesian model selection with competing sets of inequality and/or equality constraints on multiple partial correlations.

Structure Learning (Bayes Factor)
---------------------------------

### Visualizing Scientific Expectations

For Bayesian hypothesis testing in particular, it is important to *fully* understand the prior distribution ℋ<sub>*u*</sub>. This is because it captures the predicted effect size, and it is used to compute the Bayes factor. This stands in contrast to the above estimation based methods, where *E* is determined with respect to (only) the posterior distribution. To this end, **BGGM** includes functions to visualize the prior distribution or prediction--i.e.,

``` r
# define (potentially) hypothesized standard deviations for rho_ij
rho_sd <- c(0.1, 0.25, 0.5)

# plot
plt_4A <- hypothesis_plot(rho_sd = rho_sd) +
                theme(panel.grid.major = element_blank()) +
                ylab("Density")        

plt_4A
```

<img src="man/figures/README-unnamed-chunk-21-1.png" width="60%" style="display: block; margin: auto;" />

`rho_sd = c(.)` defines a couple (prospective) values for the hypothesized standard deviation of the partial correlations. Further, **BGGM** allows for testing edge differences with the Bayes factor. Accordingly, the *implied* prior distribution for the difference can also be visualized. This is implemented with:

``` r
# define (potentially) hypothesized standard deviations for rho_ij
rho_sd <- c(0.1, 0.25, 0.5)

# plot
plt_4B <- hypothesis_plot(rho_sd = rho_sd,
                          difference = TRUE) +
                theme(panel.grid.major = element_blank()) +
                ylab("Density")        

plt_4B
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="60%" style="display: block; margin: auto;" />

### Two-Sided Testing

Two-sided hypothesis testing is implemented as follows. First the model is fitted with:

``` r
# p = 5
Y <- BGGM::bfi[,1:5]

# fit model
fit_bf <- explore(Y, prior_sd = 0.5, 
                  iter = 5000, 
                  cores = 2)
summary(fit_bf)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Hypothesis Testing (Exploratory) 
#> Posterior Samples: 5000 
#> Observations (n): 2709 
#> Variables (p): 5 
#> Edges: 10 
#> Delta: 3 
#> --- 
#> Call: 
#> explore.default(X = Y, prior_sd = 0.5, iter = 5000, cores = 2)
#> --- 
#> Date: Thu May 16 15:51:31 2019
```

Note `summary(.)`, or alternatively `print(.)`, provides information about the fitted model, including that hypothesis testing (exploratory) was used, the number of observations (*n*) and variables (*p*), and the number of edges. Delta (*δ*) is the hyperparameter of the matrix−*F* distribution. A value of 3 corresponds to `prior_sd = 0.5`. This output parallels the estimation based methods. Importantly, all fitted objects include specific (what method was used) and general information (e.g., *n* and *p*) when printed.

The graph is then selected with:

``` r
E <- select(fit_bf, 
            BF_cut = 3, 
            alternative = "two.sided")

# summary 
summary(E, hyp = "H1")
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Hypothesis Testing 
#> Alternative: two.sided 
#> Bayes Factor: 3 
#> Connectivity: 80 % 
#> --- 
#> Call:
#> select.explore(x = fit_bf, BF_cut = 3, alternative = "two.sided")
#> --- 
#> Hypothesis: 
#> H1: rho != 0 
#> --- 
#> Partial Correlations 
#>  
#>            1          2          3         4         5
#> 1  0.0000000 -0.2401165 -0.1075406 0.0000000 0.0000000
#> 2 -0.2401165  0.0000000  0.2861227 0.1647403 0.1568780
#> 3 -0.1075406  0.2861227  0.0000000 0.1782371 0.3584115
#> 4  0.0000000  0.1647403  0.1782371 0.0000000 0.1210166
#> 5  0.0000000  0.1568780  0.3584115 0.1210166 0.0000000
#> --- 
#>  
#> Adjancency (non-zero) 
#>  
#>   1 2 3 4 5
#> 1 0 1 1 0 0
#> 2 1 0 1 1 1
#> 3 1 1 0 1 1
#> 4 0 1 1 0 1
#> 5 0 1 1 1 0
#> ---
```

It is also possible to change `hyp = "H1"` to `hyp = "H0`, which will print the adjacency matrix for conditionally independent relations. Further, it is possible to summarize *E* as follows:

``` r
summary(E, summarize = T, log = T, digits = 2)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Hypothesis Testing 
#> Alternative: two.sided 
#> Bayes Factor: 3 
#> --- 
#> Call:
#> select.explore(x = fit_bf, BF_cut = 3, alternative = "two.sided")
#> --- 
#> Hypotheses: 
#> H0: rho = 0
#> H1: rho != 0 
#> --- 
#> Estimates: 
#>  
#>   edge post_mean post_sd BF 10
#>  1--2   -0.2401   0.018  75.9
#>  1--3   -0.1075   0.019  12.3
#>  2--3    0.2861   0.018 111.4
#>  1--4   -0.0074   0.019  -3.5
#>  2--4    0.1647   0.019  34.1
#>  3--4    0.1782   0.019  40.3
#>  1--5   -0.0086   0.019  -3.5
#>  2--5    0.1569   0.018  32.7
#>  3--5    0.3584   0.017 191.5
#>  4--5    0.1210   0.019  16.2
#> --- 
#> note: BF_10 is evidence in favor of H1
```

The option `log = TRUE` controls the scale of the Bayes factor. Note that the same plotting options available for these methods. They have the same implementation as the estimation based methods, and thus the `alternative = two.sided` plot is not shown here. It is also possible to plot ℋ<sub>*u*</sub> and a selected posterior distribution. This visualizes how the Bayes factor is computed--i.e.,

``` r
plt_4C <- hypothesis_plot(fit = fit_bf, 
                edge = "1--4", size = 3) +
        theme(panel.grid.minor = element_blank(),
              legend.position = "top") +
  ylab("Density") +
  xlab("Distributions")


plt_4C
```

<img src="man/figures/README-unnamed-chunk-26-1.png" width="60%" style="display: block; margin: auto;" />

Here it can be seen that the Bayes factor is the ratio of density evaluated at 0. In this case, there is evidence for ℋ<sub>0</sub> (*B**F*<sub>01</sub> ≈ 33). This plotting option may be useful for understanding (and describing) the Bayes factor approach for selecting the graph (or more generally the testing strategy for partial correlations).

### One-Sided Testing

One-sided hypothesis testing is implemented with:

``` r
# p =  10
Y <- BGGM::bfi[,1:10]

# sample from posterior
fit_bf <- explore(Y, prior_sd = 0.5, 
                  iter = 5000, 
                  cores = 2)

# rho > 0
E_pos <- select(fit_bf, 
            BF_cut = 3, 
            alternative = "greater")

# positive plot
plt_pos <- plot(E_pos, type = "network") 
plt_pos <- plt_pos$plot_nonzero + 
           ggtitle(expression(atop(H[0]: rho[i][j]*" = "*0, H[1]: rho[i][j]*" > "*0))) 

plt_pos
```

<img src="man/figures/README-unnamed-chunk-27-1.png" width="60%" style="display: block; margin: auto;" />

ℋ<sub>1</sub> : *ρ*<sub>*i**j*</sub> &lt; 0 can be tested by changing `alternative = greater` to `alternative = less`.

### Exhaustive Hypothesis Testing

A defining feature of Bayesian hypothesis testing is the ability to assess which theoretical model best predicts the data at hand. However, in the absence of guiding theory, it is likely that a more exploratory approach is warranted. **BGGM** thus includes an exhaustive approach--i.e.,ℋ<sub>0</sub> : *ρ*<sub>*i**j*</sub> = 0 vs. ℋ<sub>1</sub> : *ρ*<sub>*i**j*</sub> &gt; 0 vs. ℋ<sub>2</sub> : *ρ*<sub>*i**j*</sub> &lt; 0.

which covers the entire parameter space. Further details can be found in Williams and Mulder (2019). The exhaustive approach is implemented with:

``` r
# p = 10
Y <- BGGM::bfi[,1:10]

# sample from posterior
fit_bf <- explore(Y, prior_sd = 0.5, 
                  iter = 5000, 
                  cores = 2)
# select the graph
E <- select(fit_bf, 
            hyp_prob = 0.90, 
            alternative = "exhaustive")

# first 5 rows
head(E, summarize = T, nrow = 5)
#> BGGM: Bayesian Gaussian Graphical Models 
#> --- 
#> Type: Hypothesis Testing 
#> Alternative: exhaustive 
#> Posterior probability: 0.9 
#> --- 
#> Call:
#> select.explore(x = fit_bf, alternative = "exhaustive", hyp_prob = 0.9)
#> --- 
#> Hypotheses: 
#> H0: rho = 0
#> H1: rho > 0
#> H2: rho < 0 
#> --- 
#> Estimates: 
#>  
#>   edge   post_mean    post_sd p(H0|Y) p(H1|Y) p(H2|Y)
#>  1--2 -0.24397944 0.01814374   0.000   0.000   1.000
#>  1--3 -0.10607145 0.01900275   0.000   0.000   1.000
#>  2--3  0.28679306 0.01827821   0.000   1.000   0.000
#>  1--4 -0.01427348 0.01927410   0.997   0.001   0.002
#>  2--4  0.16094706 0.01850594   0.000   1.000   0.000
#> ---
```

The ratio of posterior probabilities, e.g., *p*(ℋ<sub>1</sub>|**Y**)/*p*(ℋ<sub>2</sub>|**Y**), can then be used to compute a Bayes factor if desired. `hyp_prob = 0.90` is the decision rule for concluding there is evidence for a given hypothesis. This can the be plotted with:

``` r
plts <- plot(E, type = "network")

# p(H0 | Y) > 0.90
plts$plot_H0 +
  ggtitle(expression("p"*"("*H[0]*"| Y)") )
```

<img src="man/figures/README-unnamed-chunk-29-1.png" width="60%" style="display: block; margin: auto;" />

``` r

# p(H1 | Y) > 0.90
plts$plot_H1 +
  ggtitle(expression("p"*"("*H[1]*"| Y)") )
```

<img src="man/figures/README-unnamed-chunk-29-2.png" width="60%" style="display: block; margin: auto;" />

``` r

# p(H2 | Y) > 0.90
plts$plot_H2 +
  ggtitle(expression("p"*"("*H[2]*"| Y)") )
```

<img src="man/figures/README-unnamed-chunk-29-3.png" width="60%" style="display: block; margin: auto;" />

Comparing GGMs
==============

References
==========
