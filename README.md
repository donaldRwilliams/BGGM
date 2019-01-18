
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

There is a direct correspondence between the precision matrix, that is the inverse of the covariance matrix, and multiple regression. The details are provided ([here](https://donaldrwilliams.github.io/post/2019-10-1-r-markdown/)). Rather than fit a sequence of regression models (i.e., neighborhood selection), as in the R package [GGMnonreg](https://github.com/donaldRwilliams/GGMnonreg), it is possible to only estimate the precision matrix and then transform the elements to their respective regression counterparts. This approach is described in Williams (2018).

With the regression coefficients in-hand, it is then possible to compute *R*<sup>2</sup> for each node in the network. Similar approaches are sometimes used in the social-behavioral sciences. Here the GGMs are often estimated with ℓ<sub>1</sub>-regularization and the reported *R*<sup>2</sup> is a point estimate. This is problematic, because it can be misleading to note that one node has higher *R*<sup>2</sup> than another when there is not a measure of uncertainty. The present methods overcome this limitation with Bayesian *R*<sup>2</sup> that is described ([here](http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf))

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
| A1  |   0.1776733|  0.0166916|  0.1338691|  0.2217138|
| A2  |   0.3902242|  0.0163523|  0.3453360|  0.4306611|
| A3  |   0.4356522|  0.0142019|  0.3991760|  0.4752381|
| A4  |   0.2955884|  0.0184929|  0.2477679|  0.3355173|
| A5  |   0.3734328|  0.0172916|  0.3297305|  0.4150879|
| C1  |   0.2675226|  0.0163354|  0.2299629|  0.3083797|
| C2  |   0.3627923|  0.0134015|  0.3260684|  0.3943261|
| C3  |   0.2524328|  0.0160989|  0.2146843|  0.2965287|
| C4  |   0.3878328|  0.0151072|  0.3539124|  0.4290049|
| C5  |   0.3618108|  0.0175442|  0.3159724|  0.4033819|
| E1  |   0.3568666|  0.0155676|  0.3180822|  0.3982746|
| E2  |   0.4102940|  0.0176566|  0.3647664|  0.4548355|
| E3  |   0.4030540|  0.0154979|  0.3580366|  0.4371835|
| E4  |   0.4446836|  0.0152517|  0.4065254|  0.4807344|
| E5  |   0.3233709|  0.0179872|  0.2732149|  0.3675481|
| N1  |   0.5951378|  0.0086867|  0.5723302|  0.6145132|
| N2  |   0.5338299|  0.0123604|  0.4992442|  0.5654270|
| N3  |   0.4712103|  0.0119513|  0.4388482|  0.5019671|
| N4  |   0.4052627|  0.0177184|  0.3576245|  0.4447350|
| N5  |   0.3073446|  0.0181707|  0.2577591|  0.3493542|
| O1  |   0.2461845|  0.0143216|  0.2098434|  0.2831332|
| O2  |   0.1959535|  0.0154468|  0.1565887|  0.2335757|
| O3  |   0.3184158|  0.0151225|  0.2800979|  0.3568682|
| O4  |   0.1626688|  0.0146944|  0.1299435|  0.2040943|
| O5  |   0.2170027|  0.0153065|  0.1763381|  0.2584381|

Note that the R package mgm can also compute *R*<sup>2</sup>, but only provides point estimates. It would be possible to apply the bootstrap. However, then *R*<sup>2</sup> would not be conditional on the fitted model and posterior distributions but instead would capture sampling variability.

The package BGGM can also plot the results.

``` r
plot_r2 <- BGGM::predictive_plot(r2, 
                                 size = 2, 
                                 color = "blue")
plot_r2
```

![](man/figures/README-unnamed-chunk-7-1.png)
