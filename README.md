
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

With the regression coefficients in-hand, it is then possible to compute *R*<sup>2</sup> for each node in the network. Similar approaches are sometimes used in the social-behavioral sciences. Here the GGMs are often estimated with ℓ<sub>1</sub>-regularization and the reported *R*<sup>2</sup> is a point estimate. This is problematic, because it can be misleading to note that one node has higher *R*<sup>2</sup> than another when there is not a measure of uncertainty. The present methods overcome this limitation with Bayesian *R*<sup>2</sup> that is described [here](http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf)

JAGS will need to be installed to estimate this model ([link](https://sourceforge.net/projects/mcmc-jags/files/))

The following code fits the model, then selects the graph, and finally computes Bayesian *R*<sup>2</sup> conditional on the fitted model:

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
| A1  |   0.1883695|  0.0169427|  0.1483073|  0.2301100|
| A2  |   0.3908436|  0.0151503|  0.3519957|  0.4257576|
| A3  |   0.4247919|  0.0140904|  0.3902738|  0.4594377|
| A4  |   0.2960216|  0.0184537|  0.2493414|  0.3439533|
| A5  |   0.3726663|  0.0179936|  0.3349501|  0.4151132|
| C1  |   0.2676955|  0.0161162|  0.2250139|  0.3106030|
| C2  |   0.3612206|  0.0133233|  0.3281061|  0.3940112|
| C3  |   0.2534096|  0.0151985|  0.2185585|  0.2934040|
| C4  |   0.4000672|  0.0145855|  0.3593509|  0.4373497|
| C5  |   0.3612097|  0.0181928|  0.3186064|  0.4001974|
| E1  |   0.3565906|  0.0143239|  0.3210204|  0.3918290|
| E2  |   0.4338273|  0.0153173|  0.3922678|  0.4703743|
| E3  |   0.4102626|  0.0148775|  0.3696444|  0.4423100|
| E4  |   0.4556025|  0.0142341|  0.4169238|  0.4857838|
| E5  |   0.3238241|  0.0185447|  0.2761092|  0.3686461|
| N1  |   0.5946657|  0.0086407|  0.5706682|  0.6132380|
| N2  |   0.5337332|  0.0131920|  0.5012615|  0.5659472|
| N3  |   0.4699821|  0.0118583|  0.4395410|  0.4991964|
| N4  |   0.4256822|  0.0169256|  0.3861438|  0.4671784|
| N5  |   0.3071237|  0.0193187|  0.2551863|  0.3503473|
| O1  |   0.2437566|  0.0157097|  0.2075506|  0.2789982|
| O2  |   0.1962032|  0.0150615|  0.1563283|  0.2346939|
| O3  |   0.3027874|  0.0153374|  0.2655831|  0.3369440|
| O4  |   0.1621969|  0.0146192|  0.1181268|  0.1991478|
| O5  |   0.2168877|  0.0147817|  0.1827811|  0.2544819|

Note that the R package mgm can also compute *R*<sup>2</sup>, but only provides point estimates. It would be possible to apply the bootstrap. However, then *R*<sup>2</sup> would not be conditional on the fitted model and posterior distributions but instead would capture sampling variability.

The package BGGM can also plot the results:

``` r
plot_r2 <- BGGM::predictive_plot(r2, 
                                 size = 2, 
                                 color = "blue")
plot_r2
```

![](man/figures/README-unnamed-chunk-7-1.png)
