
<!-- README.md is generated from README.Rmd. Please edit that file -->
BGGM
====

This package is described in Williams and Mulder (2019) and Williams (2018). The methods are separated into two Bayesian approaches for inference: hypothesis testing and estimation. The former is described in Williams and Mulder (2018a), and allows for testing for the presence of edges with the Bayes factor. One-sided hypothesis testing is also possible. These methods can also provide evidence for the null hypothesis. There are extensions for confirmatory hypothesis testing in GGMs, that can include inequality or equality constraints on the partial correlations.

The estimation based methods are described in Williams (2018). The methods offer advantages compared to classical methods, in that a measure of uncertainty is provided for all parameters. For example, each node has a distribution for the variance explained (i.e., Bayesian *R*<sup>2</sup>). Measures of out-of-sample performance are also available, which also have a measure of uncertainty. The model is selected with credible interval exclusion of zero.

Williams, D. R. (2018, September 20). Bayesian Inference for Gaussian Graphical Models: Structure Learning, Explanation, and Prediction. ([pre-print](https://doi.org/10.31234/osf.io/x8dpr))

Williams, D. R., & Mulder, J. (2019, January 14). Bayesian Hypothesis Testing for Gaussian Graphical Models:Conditional Independence and Order Constraints. ([pre-print](https://doi.org/10.31234/osf.io/ypxd8))

Installation
------------

You can install BGGM from git hub with:

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/BGGM")
```

Example
-------

### Exploratory Hypothesis Testing

These methods allow for gaining evidence for both conditional dependence (*ρ* ≠ 0) and independence (*ρ* = 0). Note that GGMs are often thought to characterize conditional independence structures, but evidence for the null hypothesis of no effect is not (typically) assessed. The following uses the Bayes factor to assess evidence for the null vs. alternative hypothesis.

``` r
library(BGGM)
dat <- BGGM::bfi

# fit the exploratory approach
fit <- BGGM::bayes_explore(X = dat)

# select the network (threshold of 3)
select_graph <- BGGM::explore_select(fit, 
                                     threshold = 3, 
                                     type = "two_sided")
qgraph::qgraph(select_graph$partial_mat)
```

![](man/figures/README-example-1.png)

Some of the methods rely on sampling, so we found it most convenient to select the model after fitting-thus changing the threshold does not require refitting the model. This particular method does not require sampling from the prior or posterior distributions, but does rely on assuming normality. In the paper, Williams and Mulder (2019), it was shown these approximations performed well: the Bayes factor was consistent for model selection and invariant to the scale of the data. Sampling is possible for those not happy-and have some patience-with the normal approximation.

``` r
# select the network (threshold of 10)
select_graph <- BGGM::explore_select(fit, 
                                     threshold = 10, 
                                     type = "two_sided")
qgraph::qgraph(select_graph$partial_mat)
```

![](man/figures/README-unnamed-chunk-2-1.png)

It is likely that there is an expected direction. That is, maybe it does not make theoretical sense to have negative effects. At this time it is only possible to assume all relations are in the same direction, but this will be changed soon. One-sided hypothesis testing can be performed as follows:

``` r
# select the network (threshold of 10; one-sided)
select_graph <- BGGM::explore_select(fit, 
                                     threshold = 10, 
                                     type = "greater_than")
qgraph::qgraph(select_graph$partial_mat)
```

![](man/figures/README-unnamed-chunk-3-1.png)

Note that all the effects are now positive (i.e., the color green).

To date, the conditional independence structure of personality has not been directly assessed. Let us examine for which relations there is evidence for the null hypothesis.

``` r
# select the network (threshold of 3; two-sided)
select_graph <- BGGM::explore_select(fit, 
                                     threshold = 3, 
                                     type = "two_sided")
qgraph::qgraph(select_graph$BF_null_adj, layout = "circle")
```

![](man/figures/README-unnamed-chunk-4-1.png)

We are currently thinking of ways to plot the conditional independence structure (open to suggestions), but for now are using only the adjacency matrix.

Finally, for those interested in the substantive aspect of these networks, please see the **psych** package for the variable descriptions.

### Bayesian R-squared

There is a direct corresposnade between the precision matrix, that is the inverse of the covariance matrix, and multiple regression. The details are provided ([here](https://donaldrwilliams.github.io/post/2019-10-1-r-markdown/)). Rather than fit a sequence of regression models (i.e., neighborhood selection), as in the R packaage [GGMnonreg](https://github.com/donaldRwilliams/GGMnonreg), it is possible to only estimate the precison matrix and then transform the elements to their respective regression counterparts. This approach is described in Williams (2018).

With the regression coefficients in-hand, it is then possible to compute *R*<sup>2</sup> for each node in the network. Similar approaches are sometimes used in the social-behavioral sciences. Here the GGMs are often estimated with ℓ<sub>1</sub>-regularization and the reported *R*<sup>2</sup> is a point estimate. This is problematic, because it can be misleading to note that one node has higher *R*<sup>2</sup> than another when there is not a measure of uncertainty. The present methods overcome this limitation with computing Bayesian *R*<sup>2</sup> that is described ([here](http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf))

JAGS will need to be installed to estimate this model ([link](https://sourceforge.net/projects/mcmc-jags/files/))

The following code fits the model, then selects the graph, and finally computes Bayesian *R*<sup>2</sup> conditional on the fitted model.

``` r
# fit the model
fit <- BGGM::bayes_estimate(dat)

# select the graph
select <- BGGM::estimate_select(fit, ci_width = .99)

# Bayesian R2
r2 <- BGGM::bayes_R2(fit = fit, ci_width = .99, 
               selected = select$adjacency_mat, 
               samples = 500)
```

|     |  post\_mean|   post\_sd|       0.5%|      99.5%|
|-----|-----------:|----------:|----------:|----------:|
| A1  |   0.1762433|  0.0153795|  0.1343879|  0.2167329|
| A2  |   0.3787972|  0.0154766|  0.3384944|  0.4218465|
| A3  |   0.4243752|  0.0146088|  0.3837946|  0.4602597|
| A4  |   0.2930627|  0.0186739|  0.2403309|  0.3371595|
| A5  |   0.3723416|  0.0172815|  0.3277725|  0.4110295|
| C1  |   0.2674264|  0.0168821|  0.2237354|  0.3075276|
| C2  |   0.3626886|  0.0141462|  0.3241912|  0.3961935|
| C3  |   0.2533925|  0.0163347|  0.2139193|  0.2940771|
| C4  |   0.3875414|  0.0151338|  0.3495201|  0.4262289|
| C5  |   0.3619965|  0.0170042|  0.3121257|  0.4014649|
| E1  |   0.3420622|  0.0153976|  0.3031188|  0.3802848|
| E2  |   0.4112181|  0.0163961|  0.3691747|  0.4577512|
| E3  |   0.4028405|  0.0150143|  0.3636495|  0.4405028|
| E4  |   0.4439218|  0.0151850|  0.4115520|  0.4843417|
| E5  |   0.3250965|  0.0189092|  0.2803104|  0.3742529|
| N1  |   0.5956285|  0.0087006|  0.5679209|  0.6156081|
| N2  |   0.5325689|  0.0125474|  0.5013309|  0.5611329|
| N3  |   0.4709537|  0.0131091|  0.4371265|  0.4999997|
| N4  |   0.4055708|  0.0176145|  0.3576690|  0.4531165|
| N5  |   0.3076864|  0.0177331|  0.2571552|  0.3533014|
| O1  |   0.2454228|  0.0157277|  0.2112855|  0.2815676|
| O2  |   0.1969441|  0.0154500|  0.1625608|  0.2358583|
| O3  |   0.3027144|  0.0149969|  0.2647663|  0.3398599|
| O4  |   0.1615788|  0.0148849|  0.1254582|  0.2015237|
| O5  |   0.2168135|  0.0141829|  0.1823249|  0.2552581|

Note that the R package mgm also provides *R*<sup>2</sup>, but only a point estimate. It would be possible to apply the boostrap, but then *R*<sup>2</sup> would not conditional on this fitted model but instead would capture variability in the selected model.

The package BGGM can also plot the results.

``` r
plot_r2 <- BGGM::predictive_plot(r2, 
                                 size = 2, 
                                 color = "blue")
plot_r2
```

![](man/figures/README-unnamed-chunk-7-1.png)
