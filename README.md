
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
| A1  |   0.1843384|  0.0170868|  0.1407121|  0.2268644|
| A2  |   0.3793222|  0.0154584|  0.3401976|  0.4158084|
| A3  |   0.4251098|  0.0140937|  0.3892133|  0.4593636|
| A4  |   0.2939574|  0.0180705|  0.2474599|  0.3368746|
| A5  |   0.3733821|  0.0168506|  0.3276465|  0.4175422|
| C1  |   0.2808090|  0.0165221|  0.2372187|  0.3220607|
| C2  |   0.3624453|  0.0138960|  0.3243861|  0.3937191|
| C3  |   0.2531259|  0.0154444|  0.2131810|  0.2914723|
| C4  |   0.3920481|  0.0147570|  0.3552214|  0.4243922|
| C5  |   0.3611309|  0.0185880|  0.3146360|  0.4033821|
| E1  |   0.3418326|  0.0153472|  0.3054234|  0.3766458|
| E2  |   0.4112106|  0.0180507|  0.3629777|  0.4543671|
| E3  |   0.4124500|  0.0149382|  0.3727128|  0.4452719|
| E4  |   0.4441843|  0.0159870|  0.3997374|  0.4816921|
| E5  |   0.3237211|  0.0193466|  0.2759184|  0.3736665|
| N1  |   0.5950395|  0.0085594|  0.5717469|  0.6165439|
| N2  |   0.5332315|  0.0134098|  0.5008989|  0.5649439|
| N3  |   0.4713478|  0.0115029|  0.4415002|  0.4987279|
| N4  |   0.4046368|  0.0184112|  0.3605165|  0.4511908|
| N5  |   0.3091847|  0.0183291|  0.2555911|  0.3519177|
| O1  |   0.2450947|  0.0160598|  0.2068505|  0.2847126|
| O2  |   0.1967598|  0.0146717|  0.1606702|  0.2344459|
| O3  |   0.3159447|  0.0146435|  0.2800873|  0.3493020|
| O4  |   0.1619399|  0.0150331|  0.1174297|  0.2018222|
| O5  |   0.2167016|  0.0141681|  0.1822034|  0.2524083|

Note that the R package mgm also provides *R*<sup>2</sup>, but only provides point estimate. It would be possible to apply the bootstrap. However, then *R*<sup>2</sup> would not be conditional on this fitted model and posterior distributions but instead would capture sampling variability.

The package BGGM can also plot the results.

``` r
plot_r2 <- BGGM::predictive_plot(r2, 
                                 size = 2, 
                                 color = "blue")
plot_r2
```

![](man/figures/README-unnamed-chunk-7-1.png)
