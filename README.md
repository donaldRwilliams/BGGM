
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- ![](https://github.com/donaldRwilliams/BGGM/blob/master/man/figures/logo2.jpg) -->
This package includes methods introduced in Williams, Rast, Pericchi, and Mulder (2019), Williams and Mulder (2019), and Williams (2018). The package is built around two Bayesian approaches for inference: estimation and hypothesis testing.

The estimation based methods are described in Williams (2018). They offer advantages compared to classical methods, in that a measure of uncertainty is provided for all parameters. For example, each node has a distribution for variance explained (i.e., Bayesian *R*<sup>2</sup>). Measures of out-of-sample prediction error are available, which also have a measure of uncertainty. The model is selected with credible interval exclusion of zero or a region of practical equivalence. This allows for computing the posterior probability for an assumed *null* region--i.e., conditional independence. It is also possible to compare partial correlations.

The hypothesis testing based methods are described in Williams and Mulder (2019), and allow for testing edges (i.e., partial correlations) with the Bayes factor. One-sided hypothesis testing is also possible. These methods provide (relative) evidence for the null hypothesis. There are extensions for **confirmatory hypothesis testing** in GGMs--e.g., inequality or equality constraints on the partial correlations. This allows for comparing theoretically informed models with Bayesian model selection.

Further, it is possible to assess differences as well as similarities (i.e., the null hypothesis) between GGMs. These method were introduced in Williams, Rast, Pericchi, and Mulder (2019). Graphs are compared either with the posterior predictive distribution or Bayesian model selection. The latter allows for testing hypothesized changes in graphical structures between, for example, control and treatment groups. The posterior predictive approach is based on KL-divergence. It allows for testing the assumed (null) model of group equality for the entire graph or specific variables. These methods can be used to compare any number of GGMs.

Williams, D. R. (2018). Bayesian Inference for Gaussian Graphical Models: Structure Learning, Explanation, and Prediction. ([pre-print](https://doi.org/10.31234/osf.io/x8dpr))

Williams, D. R., & Mulder, J. (2019). Bayesian Hypothesis Testing for Gaussian Graphical Models:Conditional Independence and Order Constraints. ([pre-print](https://doi.org/10.31234/osf.io/ypxd8))

Williams, D. R., Rast, P., Pericchi, L. R., & Mulder, J. (2019). Comparing Gaussian Graphical Models with the Posterior Predictive Distribution and Bayesian Model Selection. ([pre-print](https://psyarxiv.com/yt386/))

Williams, D. R., & Mulder, J. (2019). BGGM: A R Package for Bayesian Gaussian Graphical Models. ([pre-print](https://psyarxiv.com/3b5hf/))

Outline
-------

This README is organized as follows:

-   [Installation](#installation)

-   [Estimation](#estimation)

    -   [Structure Learning (Estimation)](#structure-learning-estimation)

    -   [Edge Differences (Estimation)](#edge-differences)

    -   [Prediction](#prediction)

        -   [Bayesian *R*<sup>2</sup>](#bayesian-r2)

        -   [Leave-One-Out Cross-Validation](#leave-one-out-cross-validation)

-   [Hypothesis Testing](#hypothesis-testing)

    -   [Structure Learning (Bayes Factor)](#structure-learning-bayes-factor)

        -   [Visualizing Scientific Expectations](#visualizing-scientific-expectations)

        -   [Two-Sided Hypothesis Testing](#two-sided-testing)

        -   [One-Sided Hypothesis Testing](#one-sided-testing)

        -   [Exhaustive Hypothesis Testing](#exhaustive-hypothesis-testing)

    -   [Edge Differences (Hypothesis Testing)](#edge-differences-hypothesis-testing)

    -   [Confirmatory Hypothesis Testing](#confirmatory-hypothesis-testing)

        -   [Order Constraints](#order-constraints)

        -   [Equality Constraints](#equality-constraints)

-   [Comparing GGMs](#comparing-ggms)

    -   [Posterior Predictive KL Divergence](#posterior-predictive-kl-divergence)

    -   [Bayesian Hypothesis Testing](#bayesian-hypothesis-testing)

        -   [Exploratory Hypothesis Testing](#exploratory-hypothesis-testing)

        -   [Group Confirmatory Hypothesis Testing](#group-confirmatory-hypothesis-testing)

Installation
============

Note in preparation for submitting to CRAN, we ran into some issues with the submission guidelines. Hence, we are currently going through the functions one by one and finalizing them for a stable release.

You can install BGGM from git hub with.

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/BGGM")
```

<!-- # Estimation -->
<!-- The following methods are described in Williams (2018). They build upon two basic ideas. First, the Wishart distribution is a conjugate prior distribution for the precision matrix (inverse of the covariance matrix). This provides an analytic solution for selecting the graph, and allows for conveniently drawing posterior samples. Second, there is an exact relationship between estimating the precision matrix directly and with multiple regression. Here the individual elements, from joint posterior distribution of the precision matrix, can be converted to their respective regression counterparts. This allows for assessing nodewise (for each variable in the model) predictability.  -->
<!-- ##  Structure Learning (Estimation) -->
<!-- By structure learning we are referring to selecting the graph (i.e., the edge set $E$), which consists of those edges determined to be non-zero.  For demonstrative purposes, we consider a relatively small number of variables ($p= 5$).  -->
<!-- The package **BGGM** offers a convenient analytic solution for estimating GGMs. It is implemented with: -->
<!-- ```{r eval=TRUE, message=F, warning=F} -->
<!-- library(BGGM) -->
<!-- library(ggplot2) -->
<!-- library(ggraph) -->
<!-- library(foreach) -->
<!-- # p = 5 -->
<!-- Y <- BGGM::bfi[,1:5] -->
<!-- # analytic solution -->
<!-- fit_analytic <- estimate(Y, analytic = T) -->
<!-- # summary -->
<!-- summary(fit_analytic) -->
<!-- ``` -->
<!-- Note `summary(.)` provides information about the fitted model, including that the analytic solution was used, the number of observations ($n$) and variables ($p$), and the number of edges.  -->
<!-- The edge set is then selected with: -->
<!-- ```{r } -->
<!-- # select the graph (edge set E) -->
<!-- E <- select(fit_analytic, ci_width = 0.95) -->
<!-- # summary of E -->
<!-- summary(E) -->
<!-- ``` -->
<!-- The analytic solution works directly with the precision matrix, and thus, there is not an option to summarize the posterior distributions. This is because the non-standardized elements are in the opposite direction ($\pm$) of the partial correlations, which in our experience, can lead to confusion. To summarize  the partial correlations change `analytic = T` to `analytic = F`: -->
<!-- ```{r} -->
<!-- # sample from posterior -->
<!-- fit_sampling <- estimate(Y, analytic = F) -->
<!-- # select the graph -->
<!-- E <- select(fit_sampling, ci_width = 0.95) -->
<!-- # summarize partial correlations -->
<!-- summary(E, summarize = T, digits = 2) -->
<!-- ``` -->
<!-- Note that `edge` corresponds to that particular entry in the partial correlation matrix--i.e., `1--2` is the relation between the first and second variables, respectively. -->
<!-- **BGGM** provides several options for plotting, with each implemented as a S3 generic. For example, the partial correlations can be plotted with: -->
<!-- ```{r message=F, fig.height=3, fig.width=8} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # sampling required -->
<!-- fit_sampling <- estimate(Y, analytic = F) -->
<!-- # plot -->
<!-- plot_1A <- plot(fit_sampling,  -->
<!--                 ci_width = 0.95,  -->
<!--                 width = 0.1,   -->
<!--                 size = 2) + -->
<!--             coord_cartesian() + -->
<!--             theme(axis.text.x = element_text(angle = 90)) -->
<!-- plot_1A -->
<!-- ``` -->
<!-- This example nicely demonstrates how the `plot` objects can be further customized with **ggplot2**. There are two options for visualizing the selected graph. The heatmap plot is generated with: -->
<!-- ```{r, message=F} -->
<!-- # select the graph -->
<!-- E <- select(fit_sampling, ci_width = 0.95) -->
<!-- # heatmap plot -->
<!-- plot_1B <- plot(E,  -->
<!--                 type = "heatmap",  -->
<!--                 lower_tri = TRUE) + -->
<!--            ggtitle("Heatmap Plot") +  -->
<!--            theme(plot.title = element_text(size = 15)) -->
<!-- ``` -->
<!-- Here `lower_tri = TRUE` controls which partial correlations are plotted.  In this case, only the lower triangular elements are included in the plot. This can be changed with `lower_tri = FALSE`. -->
<!-- On the other hand, a “network” plot can be obtained with: -->
<!-- ```{r, warning=F,  fig.width=10, fig.height=4} -->
<!-- # network plot -->
<!-- plot_1C <- plot(E, type = "network", -->
<!--                 layout ='circle', -->
<!--                 node_outer = 8, -->
<!--                 node_inner = 7, -->
<!--                 node_text_size = 4) + -->
<!--            ggtitle("Network Plot") + -->
<!--            theme(plot.title = element_text(size = 15)) -->
<!-- cowplot::plot_grid(plot_1B, plot_1C) -->
<!-- ``` -->
<!-- A  key  feature  of **BGGM** is  extending  inference  beyond  identifying  non-zero  partial correlations.  The region of practical equivalence can be used for this purpose, as it allows for determining which relations are practically zero. In this case, we follow Cohen’s guidelines, wherein 0.1 is considered a small effect.This is implemented with: -->
<!-- ```{r, eval=T} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # sample from posterior -->
<!-- fit_sample <- estimate(Y, iter = 5000, analytic = F) -->
<!-- # select the graph -->
<!-- E <- select(fit_sample, rope = 0.1, prob = 0.95) -->
<!-- # summary for first 10 rows -->
<!-- head(E, nrow = 10, summarize = T, digits = 2) -->
<!-- ``` -->
<!-- The argument `prob = 0.95` requires that 95 % of the posterior density be in or out of the rope to be considered practically equivalent or different from zero.  With this decision rule, as seen with `head(.)`, edges `1--4` and `1--5` are practically equivalent to zero.  This inference is made possible with **BGGM**. -->
<!-- In this case, `plot(.)` returns two objects:  (1) the selected edges; (2) those for which there is support for the null values.  This is implemented with: -->
<!-- ```{r warning = F , out.height= "250 %"} -->
<!-- # network plot -->
<!-- plts <- plot(E, type = "network", -->
<!--              layout ='circle', -->
<!--              node_outer = 10, -->
<!--              node_inner = 9, -->
<!--              node_text_size = 6) -->
<!-- # practically non-zero -->
<!-- plot_1D <- plts$plot_nonzero + -->
<!--              ggtitle("Practically Non-zero") + -->
<!--              theme(plot.title = element_text(size = 15)) -->
<!-- plot_1D -->
<!-- ``` -->
<!-- ```{r, warning=F, out.height= "250 %"} -->
<!-- # practically zero -->
<!-- plot_1E <- plts$plot_zero + -->
<!--               ggtitle("Practically Zero") + -->
<!--               theme(plot.title = element_text(size = 15)) -->
<!-- plot_1E -->
<!-- ``` -->
<!-- We emphasize that GGMs are often thought to capture conditionally *independent* relations--i.e., evidence for the null hypothesis of no effect, conditional on the other variables in the model. However, the dominant approach assesses conditional *dependence* ($\rho_{ij} \neq 0$), and then sets relations to zero otherwise. **BGGM** can explicitly answer the question of conditional independence. -->
<!-- ## Edge Differences -->
<!-- Differences between partial correlations are often tested in GGMs; for example, with a classical (i.e., frequentist) approach that is implemented in **bootnet**.  One contribution of **BGGM** is providing Bayesian analogs for commonly used methods, as well as extensions to those methods.  In this case, we can use posterior probabilities to determine which edges are practically equivalent.  This is implemented with: -->
<!-- ```{r, message=F} -->
<!-- edge_difference <- edge_compare(fit_sample,  -->
<!--                                 contrast =  list("1--5 - 1--3",  -->
<!--                                                  "1--2 - 1--6",  -->
<!--                                                  "1--4 - 1--7",  -->
<!--                                                  "1--5 - 1--10",  -->
<!--                                                  "1--2 - 1--9"),  -->
<!--                                 ci_width = 0.95, -->
<!--                                 rope = 0.1) -->
<!-- head(edge_difference, nrow = 4) -->
<!-- ``` -->
<!-- This  output  includes  the  posterior  mean  and  standard  deviation  for  each  difference. Further, `pr_in` is  the  proportion  of  samples  between ($\pm$) 0.1.   This  can  be interpreted as the posterior probability of practical equivalence, which has been defined with the argument `rope = 0.1`.  Note that care must be taken when specifying the contrasts, as an error will arise if they are not in the proper format. The object `edge_difference` can then be plotted, but this is omitted to save space. -->
<!-- ## Prediction -->
<!-- The following is based on the correspondence between the elements of the precision matrix and multiple regression. In the context of GGMs, using regression to select edges is referred to as “neighborhood” selection.   On the other hand, the method described in Williams (2018) works directly with either the posterior distribution for the precision matrix or the maximum a posteriori estimates. These are then converted to the corresponding regression coefficients and residual variances. It follows that **BGGM** can also be used for the purpose of multiple regression–i.e., -->
<!-- ```{r} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # sample posterior -->
<!-- fit <- estimate(Y, iter = 5000) -->
<!-- # precision to regression -->
<!-- coefficients(fit, node = 1, ci_width = 0.95) -->
<!-- ``` -->
<!-- Here `node = 1` indicates which node is summarized.  This correspondence allows for computing measures of prediction error (or accuracy), including Bayesian $R^2$ and Bayesian leave-one-out cross-validation,  each  of  which  has  a  measure  of  uncertainty. Furthermore,  when  a computationally convenient option is desirable, **BGGM** includes an analytic expression for prediction error. This is also known as the predicted residual sums of squares (PRESS). -->
<!-- ### Bayesian $R^2$ -->
<!-- In-sample Bayesian $R^2$ is implemented with: -->
<!-- ```{r} -->
<!-- # training data -->
<!-- Y_train <- BGGM::bfi[1:100,1:10] -->
<!-- # fit to training data -->
<!-- fit_train <- estimate(Y_train, iter = 5000) -->
<!-- # compute Bayes R2 -->
<!-- train_R2 <- predict(fit_train, -->
<!--                     ci_width = 0.90, -->
<!--                     samples = 1000, -->
<!--                     measure = "R2") -->
<!-- # summary for first 2 rows -->
<!-- head(train_R2, nrow = 2) -->
<!-- ``` -->
<!-- Here `ci_width = 0.90` indicates the decision rule for setting coefficients to zero, and by default, 95 % intervals are used in the summary output. Similarly, out-of-sample Bayesian $R^2$ is computed with: -->
<!-- ```{r} -->
<!-- # test data -->
<!-- Y_test <-  BGGM::bfi[101:2000,1:10] -->
<!-- # predict test data -->
<!-- test_R2 <- predict(fit_train, ci_width = 0.90, -->
<!--                    test_data = Y_test, -->
<!--                    samples = 1000, measure = "R2") -->
<!-- ``` -->
<!-- The work flow is completed by visualizing Bayesian $R^2$ for each node–i.e., -->
<!-- ```{r, out.height= "300 %"} -->
<!-- # prior training and test error in the same plot -->
<!-- plt_3A <- plot(x1 = train_R2, x2 =  test_R2, order = "test") -->
<!-- plt_3A -->
<!-- ``` -->
<!-- Here the nodes have been ordered by which has the best out-of-sample performance. It is also possible to have each in a separate plot by leaving `x2` empty. The `predict` object can be used to assess differences in predictive accuracy with compare(.). **BGGM** also includes mean squared error (`measure = "mse"`). -->
<!-- ### Leave-one-out cross-validation -->
<!-- Bayesian leave-one-out cross-validation is implemented with: -->
<!-- ```{r} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[1:1000,1:10] -->
<!-- # sample posterior -->
<!-- fit_sample <- estimate(Y, iter = 5000) -->
<!-- # Bayesian LOO -->
<!-- bayes_loo <- loocv(fit_sample) -->
<!-- # nodewise loo summary -->
<!-- summary(bayes_loo) -->
<!-- ``` -->
<!-- The results are plotted with: -->
<!-- ```{r, out.height= "300 %"} -->
<!-- # plot CV error -->
<!-- plt_3B <- plot(bayes_loo, size = 8) + -->
<!--           theme_classic() + -->
<!--           ylab("Bayesian Leave-One-Out") -->
<!-- plt_3B -->
<!-- ``` -->
<!-- Similarly,  by  setting `analytic = T`,  leave-one-out  prediction  error  can  be  computed analytically.  This is implemented with: -->
<!-- ```{r, out.height= "300 %"} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[1:1000,1:10] -->
<!-- # analytic solution -->
<!-- fit_analytic <- estimate(Y, analytic = T) -->
<!-- # analytic LOO (PRESS; based on point estimates) -->
<!-- press_loo <- loocv(fit_analytic) -->
<!-- # plot CV error -->
<!-- plt_3C <- plot(press_loo, size = 8) + -->
<!--           theme_classic() + -->
<!--           ylab("PRESS: Leave-One-Out") + -->
<!--           scale_y_continuous(expand = c(0, 0), -->
<!--           limit = c(0, 1000)) -->
<!-- plt_3C -->
<!-- ``` -->
<!-- This highlights the difference between the leave-one-out methods, in that the Bayesian version has a measure of uncertainty (although the order is the same). For both measures of predictive *error*, a lower value indicates a more predictable node (variable). -->
<!-- # Hypothesis Testing -->
<!-- The following methods were introduced in Williams and Mulder (2019). That work not only presented an exploratory approach using the Bayes factor, but it also proposed methodology for confirmatory hypothesis testing in GGMs. The latter provides an alternative to data driven model selection that is commonplace in the GGM literature, and in particular, it allows for comparing theoretical models. Further, the novel matrix$-F$ prior distribution is used for the partial correlations, which offers more flexibility than the Wishart distribution. This approach builds upon (Mulder2016), where the focus was on correlations (as opposed to *partial* correlations). In particular, **BGGM** allows for Bayesian model selection with competing sets of inequality and/or equality constraints on multiple partial correlations. -->
<!-- ## Structure Learning (Bayes Factor) -->
<!-- ### Visualizing Scientific Expectations -->
<!-- For Bayesian hypothesis testing in particular, it is important to *fully* understand the prior distribution $\mathcal{H}_u$. This is because it captures the predicted effect size, and it is used to compute the Bayes factor. This stands in contrast to the above estimation based methods, where $E$ is determined with respect to (only) the posterior distribution. To this end, **BGGM** includes functions to visualize the prior distribution or prediction--i.e., -->
<!-- ```{r,  out.width = '50%'} -->
<!-- # define (potentially) hypothesized standard deviations for rho_ij -->
<!-- rho_sd <- c(0.1, 0.25, 0.5) -->
<!-- # plot -->
<!-- plt_4A <- hypothesis_plot(rho_sd = rho_sd) + -->
<!--                 theme(panel.grid.major = element_blank()) + -->
<!--                 ylab("Density")         -->
<!-- plt_4A -->
<!-- ``` -->
<!-- `rho_sd = c(.)` defines a couple (prospective) values for the hypothesized standard deviation of the partial correlations. Further, **BGGM** allows for testing edge differences with the Bayes factor. Accordingly, the *implied* prior distribution for the difference can also be visualized. This is implemented with: -->
<!-- ```{r,  out.width = '50%'} -->
<!-- # define (potentially) hypothesized standard deviations for rho_ij -->
<!-- rho_sd <- c(0.1, 0.25, 0.5) -->
<!-- # plot -->
<!-- plt_4B <- hypothesis_plot(rho_sd = rho_sd, -->
<!--                           difference = TRUE) + -->
<!--                 theme(panel.grid.major = element_blank()) + -->
<!--                 ylab("Density")         -->
<!-- plt_4B -->
<!-- ``` -->
<!-- ### Two-Sided Testing -->
<!-- Two-sided hypothesis testing is implemented as follows. First the model is fitted with: -->
<!-- ```{r, warning=F} -->
<!-- # p = 5 -->
<!-- Y <- BGGM::bfi[,1:5] -->
<!-- # fit model -->
<!-- fit_bf <- explore(Y, prior_sd = 0.5,  -->
<!--                   iter = 5000,  -->
<!--                   cores = 2) -->
<!-- summary(fit_bf) -->
<!-- ``` -->
<!-- Note `summary(.)`, or alternatively `print(.)`, provides information about the fitted model, including that hypothesis testing (exploratory) was used, the number of observations ($n$) and variables ($p$), and the number of edges. Delta ($\delta$) is the hyperparameter of the matrix$-F$ distribution. A value of 3 corresponds to `prior_sd = 0.5`. This output parallels the estimation based methods. Importantly, all fitted objects include specific (what method was used) and general information (e.g., $n$ and $p$) when printed. -->
<!-- The graph is then selected with: -->
<!-- ```{r} -->
<!-- E <- select(fit_bf,  -->
<!--             BF_cut = 3,  -->
<!--             alternative = "two.sided") -->
<!-- # summary  -->
<!-- summary(E, hyp = "H1") -->
<!-- ``` -->
<!-- It is also possible to change `hyp = "H1"` to `hyp = "H0`, which will print the adjacency matrix for conditionally independent relations. Further, it is possible to summarize $E$ as follows: -->
<!-- ```{r} -->
<!-- summary(E, summarize = T, log = T, digits = 2) -->
<!-- ``` -->
<!-- The option `log = TRUE` controls the scale of the Bayes factor. Note that the same plotting options available for these methods. They have the same implementation as the estimation based methods, and thus the `alternative = two.sided` plot is not shown here. It is also possible to plot $\mathcal{H}_u$ and a selected posterior distribution. This visualizes how the Bayes factor is computed--i.e., -->
<!-- ```{r,  out.width = '50%'} -->
<!-- plt_4C <- hypothesis_plot(fit = fit_bf,  -->
<!--                 edge = "1--4", size = 3) + -->
<!--         theme(panel.grid.minor = element_blank(), -->
<!--               legend.position = "top") + -->
<!--   ylab("Density") + -->
<!--   xlab("Distributions") -->
<!-- plt_4C -->
<!-- ``` -->
<!-- Here it can be seen that the Bayes factor is the ratio of density evaluated at 0. In this case, there is evidence for $\mathcal{H}_0$ ($BF_{01} \approx 33$). This plotting option  may be useful for understanding (and describing) the Bayes factor approach for selecting the graph (or more generally the testing strategy for partial correlations).  -->
<!-- ### One-Sided Testing -->
<!-- One-sided hypothesis testing is implemented with: -->
<!-- ```{r, warning=F, message=F, out.width = '50%'} -->
<!-- # p =  10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # sample from posterior -->
<!-- fit_bf <- explore(Y, prior_sd = 0.5,  -->
<!--                   iter = 5000,  -->
<!--                   cores = 2) -->
<!-- # rho > 0 -->
<!-- E_pos <- select(fit_bf,  -->
<!--             BF_cut = 3,  -->
<!--             alternative = "greater") -->
<!-- # positive plot -->
<!-- plt_pos <- plot(E_pos, type = "network")  -->
<!-- plt_pos <- plt_pos$plot_nonzero +  -->
<!--            ggtitle(expression(atop(H[0]: rho[i][j]*" = "*0, H[1]: rho[i][j]*" > "*0)))  -->
<!-- plt_pos -->
<!-- ``` -->
<!-- $\mathcal{H}_1: \rho_{ij} < 0$ can be tested by changing `alternative = greater` to `alternative = less`. -->
<!-- ### Exhaustive Hypothesis Testing -->
<!-- A defining feature of Bayesian hypothesis testing is the ability to assess which theoretical model best predicts the data at hand. However, in the absence of guiding theory, it is likely that a more exploratory approach is warranted. **BGGM** thus includes an exhaustive approach--i.e.,$\mathcal{H}_0: \rho_{ij} = 0$ vs. $\mathcal{H}_1: \rho_{ij} > 0$ vs. $\mathcal{H}_2: \rho_{ij} < 0$. -->
<!-- This covers the entire parameter space. Further details can be found in Williams and Mulder (2019). The exhaustive approach is implemented with: -->
<!-- ```{r} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # sample from posterior -->
<!-- fit_bf <- explore(Y, prior_sd = 0.5,  -->
<!--                   iter = 5000,  -->
<!--                   cores = 2) -->
<!-- # select the graph -->
<!-- E <- select(fit_bf,  -->
<!--             hyp_prob = 0.90,  -->
<!--             alternative = "exhaustive") -->
<!-- # first 5 rows -->
<!-- head(E, summarize = T, nrow = 5) -->
<!-- ``` -->
<!-- The ratio of posterior probabilities, e.g., $p(\mathcal{H}_1|\textbf{Y}) / p(\mathcal{H}_2|\textbf{Y})$, can then be used to compute a Bayes factor if desired. `hyp_prob = 0.90` is the decision rule for concluding there is evidence for a given hypothesis. This can the be plotted with: -->
<!-- ```{r, warning=F, message=F, fig.width=10, fig.height=3.5} -->
<!-- # plots -->
<!-- plts <- plot(E, type = "network") -->
<!-- # combine -->
<!-- cowplot::plot_grid(plts$plot_H0, plts$plot_H1, plts$plot_H2, nrow = 1) -->
<!-- ``` -->
<!-- ## Edge Differences (Hypothesis Testing) -->
<!-- **BGGM** allows for testing edge differences with Bayesian hypothesis testing. The same options are available as for determining $E$--i.e., `alternative = "two.sided"`, `alternative = "greater"` (or less), and `alternative = "exhaustive"`. The exhaustive approach is implemented with: -->
<!-- ```{r, out.height= "350 %"} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # sample from posterior -->
<!-- fit_bf <- explore(Y, prior_sd = 0.5,  -->
<!--                   iter = 5000,  -->
<!--                   cores = 2) -->
<!-- # edge compare -->
<!-- edge_comp <- edge_compare(fit_bf,  -->
<!--                   contrast = c("1--5 - 1--3",  -->
<!--                                "1--2 - 1--6",  -->
<!--                                "1--4 - 1--7",  -->
<!--                                "1--5 - 1--10",  -->
<!--                                "1--2 - 1--9"),  -->
<!--                   alternative = "exhaustive") -->
<!-- # summary -->
<!-- summary(edge_comp) -->
<!-- ``` -->
<!-- There are two plotting options for the returned object `edge_comp`. The first is ideal if few contrasts are tested--i.e., -->
<!-- ```{r, out.height= "350 %"} -->
<!-- # plot -->
<!-- plt <- plot(edge_comp, stack = F, spread = .75) -->
<!-- # one contrast -->
<!-- plt$`1--5 - 1--3` -->
<!-- ``` -->
<!-- The object `plt` includes a separate plot for each contrast. On the other hand, in the case of many contrasts, it is also possible to visualize each with a stacked bar chart--i.e., -->
<!-- ```{r, out.height= "350 %"} -->
<!-- # plot -->
<!-- plt <- plot(edge_comp, stack = T, spread = .75) -->
<!-- # stacked -->
<!-- plt -->
<!-- ``` -->
<!-- # Confirmatory Hypothesis Testing -->
<!-- A key contribution of **BGGM** is extending hypothesis testing beyond exploratory and to confirmatory in GGMs. The former is essentially feeding -->
<!-- the data to the functions in **BGGM** and seeing what comes back. In other words, there are no specific, hypothesized models under consideration. On the other hand, the confirmatory hypothesis testing approaches allows for comparing theoretical models or (actual) predictions. The focus is thus not on $E$, but only certain edges in the network.  -->
<!-- ## Order Constraints -->
<!-- Theory may suggest, for example, that a set of partial correlations is expected to be larger than another set of partial correlations, that there is a hypothesized order of edge magnitude for a given node (variable), or some edges are expected to be positive while others are predicted to be negative. These can be tested with **BGGM**.  -->
<!-- For example, suppose we expected (predicted) the following:  -->
<!-- $\mathcal{H}_1: \rho_{1,2} > \rho_{1,3} > \rho_{1,4} > \rho_{1,5}$ -->
<!-- vs.   -->
<!-- $\mathcal{H}_c:$ Not $\mathcal{H}_1$ -->
<!-- This is implemented with: -->
<!-- ```{r} -->
<!-- # p = 10 -->
<!-- Y <- BGGM::bfi[,1:10] -->
<!-- # hypothesis -->
<!-- hypothesis <- c("1--2 > 1--3 > 1--4 > 1--5") -->
<!-- # test order -->
<!-- test_order <-  confirm(x = Y, hypothesis  = hypothesis,  -->
<!--                        prior_sd = 0.5, iter = 50000,  -->
<!--                        cores = 2) -->
<!-- # summary -->
<!-- summary(test_order) -->
<!-- ``` -->
<!-- As indicated by the posterior probabilities, the hypothesis was not supported. -->
<!-- One the other hand, we can also test the prediction that, for a specific variable (for example),  certain edges are negative while others are positive--i.e.,   -->
<!-- $\mathcal{H}_1: (\rho_{1,2}, \rho_{1,3}, \rho_{1,4}) < 0 < \rho_{1,6}$ -->
<!-- vs. -->
<!-- $\mathcal{H}_c:$ Not $\mathcal{H}_1$ -->
<!-- Here the prediction states that the relation between variable 1 and 6 is positive, whereas the others are predicted to be negative. This is implemented with:  -->
<!-- ```{r} -->
<!-- # hypothesis -->
<!-- hypothesis <- c("(1--2, 1--3, 1--4)  <  0 < (1--6)") -->
<!-- # test order -->
<!-- test_order <-  confirm(x = Y, hypothesis  = hypothesis,  -->
<!--                        prior_sd = 0.5, iter = 50000,  -->
<!--                        cores = 2) -->
<!-- # summary -->
<!-- summary(test_order) -->
<!-- ``` -->
<!-- There is strong support for $\mathcal{H}_1$. Note that `iter = 50000` ($\approx 7$ seconds), which ensures an adequate number of samples for accurately computing the Bayes factor. Although the focus is not *E*, but rather a subset of partial correlations, the object `test_order` can be plotted--i.e.,  -->
<!-- ```{r, out.width= "325 %", out.width= "325 %"} -->
<!-- plt <- hypothesis_plot(fit = test_order,  -->
<!--                 node_outer = 8,  -->
<!--                 node_inner = 7,  -->
<!--                 node_text_size = 6) -->
<!-- plt -->
<!-- ``` -->
<!-- This is meant to highlight the specific edges considered in the hypothesis test. Further, in this plot, it is clear that $(\rho_{1,2}, \rho_{1,3}, \rho_{1,4}) < 0$ and that $0 < \rho_{1,6}$. We emphasize that these *confirmatory* hypotheses should be decided upon before visualizing the network. -->
<!-- ## Equality Constraints -->
<!-- **BGGM** can also test exact equality constraints on sets of edges (i.e., partial correlations). One could predict that certain edges have the same strength (effect size) or that a set of edges are equal to zero. The latter allows for testing expectations about the underlying conditional independence structure, say, for a particular variable (node) in the model. For example, -->
<!-- $\mathcal{H}_1: (\rho_{1,2}, \rho_{1,3}, \rho_{1,4}) = 0$ -->
<!-- vs. -->
<!-- $\mathcal{H}_c:$ Not $\mathcal{H}_1$ -->
<!-- ```{r} -->
<!-- hypothesis <- "(1--2, 1--3, 1--4) = 0" -->
<!-- # test order -->
<!-- test_order <-  confirm(x = Y, hypothesis  = hypothesis,  -->
<!--                        prior_sd = 0.5, iter = 50000,  -->
<!--                        cores = 2) -->
<!-- # summary -->
<!-- summary(test_order) -->
<!-- ``` -->
<!-- # Comparing GGMs -->
<!-- **BGGM** includes two approaches for comparing any number of GGMs. They were introduced in Williams, Rast, Pericchi, and Mulder (2019). The methods can test entire structures, particular edges, and nodewise testing. One approach is based on the posterior predictive distribution, and tests the null hypothesis of group equality. Importantly, like classical methods (i.e., frequentist), this approach can only determine whether the null of group equality is rejected or retained. It cannot provide evidence for the null hypothesis. To this end, another approach is based on Bayesian model selection. This allows for assessing network invariances (i.e., the null hypothesis). -->
<!-- ## Posterior Predictive KL Divergence -->
<!-- The predictive approach is based on the distribution of future data. In particular, the posterior predictive distribution of Kullback–Leibler divergence is used to test the null hypothesis of group equality. This method can be understood as a multivariate likelihood ratio that accounts for uncertainty, but importantly, with respect to the predictive distribution. -->
<!-- ```{r} -->
<!-- # Assume null is true -->
<!-- Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1) -->
<!-- Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1) -->
<!-- Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1) -->
<!-- # predictive check -->
<!-- ggm_ppc1 <- ggm_compare_ppc(Y1, Y2, Y3,  -->
<!--                            type = "global", iter = 5000) -->
<!-- # summary -->
<!-- summary(ggm_ppc1) -->
<!-- ``` -->
<!-- In this case, the null hypothesis was not rejected. The `p_value` is the probability, or proportion, that the predictive KL divergence is larger than the observed KL divergence. For example, with `p_value = 0.02` (e.g.), this would indicate that, conditional on group equality, there is only a 2 % chance of observing that amount of error in the future. Thus this would be considered extreme at the conventional threshold ($\alpha = 0.05$). -->
<!-- The above example tested the entire network (`type = global`). It is also possible to narrow the focus to each node in the GGM. In other words, the null hypothesis of group equality is tested for each variable. This predictive check is implemented with: -->
<!-- ```{r} -->
<!-- # Assume null is true -->
<!-- Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1) -->
<!-- Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1) -->
<!-- # predictive check -->
<!-- ggm_ppc2 <- ggm_compare_ppc(Y1, Y2, type = "nodewise", iter = 5000) -->
<!-- # summary -->
<!-- summary(ggm_ppc2) -->
<!-- ``` -->
<!-- The predictive distributions can also be visualized. For example, the `global` check is implemented with: -->
<!-- ```{r, out.width= "500 %", out.width= "500 %", message=F, warning=F} -->
<!-- # global predictive check -->
<!-- plot(ggm_ppc1) + -->
<!--   theme_minimal(base_size = 14)  + -->
<!--   theme(legend.position = "none",  panel.grid.minor.x = element_blank()) -->
<!-- ``` -->
<!-- The `nodewise` predictive check is plotted with: -->
<!-- ```{r, message=F, warning=F, out.height= "500 %"} -->
<!-- #plots -->
<!-- plot(ggm_ppc2, log = T)[[1]] + -->
<!--   theme_minimal(base_size = 14)  + -->
<!--   theme(legend.position = "none",  panel.grid.minor.x = element_blank()) -->
<!-- ``` -->
<!-- ## Bayesian Hypothesis Testing -->
<!-- ### Exploratory Hypothesis Testing -->
<!-- The following method uses the Bayes factor to test whether there is (relative) evidence for edge equality in any number of GGMs. The test is for each edge or partial  -->
<!-- correlation in the GGMs. In the case of more than two groups, for example with, say, 4 groups, the method tests whether each edge is the same across all groups. The alternative ($\mathcal{H}_1$) is then the unrestricted model--i.e., not $\mathcal{H}_0$. This is implemented with  -->
<!-- ```{r} -->
<!-- # Assume null is true for all edges -->
<!-- Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16)) -->
<!-- Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16)) -->
<!-- Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16)) -->
<!-- # sample from posteriors and compute BF -->
<!-- ggm_bf <- ggm_compare_bf(Y1, Y2, Y3,  -->
<!--                          prior_sd = .5,  -->
<!--                          iter = 5000,  -->
<!--                          cores = 2) -->
<!-- # select with BF_cut = 3 -->
<!-- ggm_bf_sel <- select(ggm_bf, BF_cut = 3) -->
<!-- # summary -->
<!-- summary(ggm_bf_sel, type = "adj") -->
<!-- ``` -->
<!-- ### Group Confirmatory Hypothesis Testing -->
<!-- The confirmatory methods described above can be extended to any number of groups.  For example, a researcher can test their expectations in, say, treatment vs. control groups, where certain edges are predicted to increase. This can then be compared to a model predicting those edges stayed the same--i.e., -->
<!-- $\mathcal{H}_1: (\rho_{1,2_{g1}} > \rho_{1,2_{g2}}) ,  (\rho_{1,3_{g1}} > \rho_{1,3_{g2}}), (\rho_{2,3_{g1}} > \rho_{2,3_{g2}})$  -->
<!-- vs. -->
<!-- $\mathcal{H}_2: (\rho_{1,2_{g1}} = \rho_{1,2_{g2}}) ,  (\rho_{1,3_{g1}} = \rho_{1,3_{g2}}), (\rho_{2,3_{g1}} = \rho_{2,3_{g2}})$ -->
<!-- vs. -->
<!-- $\mathcal{H}_c:$  not $\mathcal{H}_1$ or $\mathcal{H}_2$ -->
<!-- Note that this method can be used for any number of groups, variables, and sets of competing inequality and equality restrictions. This example is implemetned with:  -->
<!-- ```{r} -->
<!-- # assume the null is true -->
<!-- Y1 <- MASS::mvrnorm(500, rep(0, 3), Sigma = diag(3)) -->
<!-- Y2 <- MASS::mvrnorm(500, rep(0, 3), Sigma = diag(3)) -->
<!-- ggms_confirm <- ggm_compare_bf(Y1, Y2, prior_sd = 0.5,  -->
<!--                                 hypothesis  = "g1_1--2 > g2_1--2,  -->
<!--                                                g1_1--3 > g2_1--3,  -->
<!--                                                g1_2--3 > g2_2--3;  -->
<!--                                                g1_1--2 = g2_1--2,  -->
<!--                                                g1_1--3 = g2_1--3,  -->
<!--                                                g1_2--3 = g2_2--3") -->
<!-- summary(ggms_confirm) -->
<!-- ``` -->
<!-- In this case, because each group was assumed to have equal edges, $\mathcal{H}_2$ was supported by the data. We again emphasize that edges can be tested in any size network (assuming that $p < n$). -->
