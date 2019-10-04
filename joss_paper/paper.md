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

Gaussian graphical models (GGM) allow for learning conditional independence structures that are encoded by partial correlations. Whereas there are several R packages for classical (i.e., frequentist) methods, there are only two that implement a Bayesian approach. These are exclusively focused on identifying the graphical structure. Thepackage **BGGM** not only contains novel Bayesian methods for this purpose, but it also includes Bayesian methodology for extending inference beyond identifying non-zero relations. **BGGM** is built around two Bayesian approaches for inference--estimation and hypothesis testing. The former focuses on the posterior distribution and includes extensions to assess predictability, as well as methodology to compare partial correlations. The latter includes methods for Bayesian hypothesis testing, in both exploratory and confirmatory contexts, with the novel matrix-$F$ prior distribution. This allows for testing order and equality constrained hypotheses, as well as a combination of both with the Bayes factor. Further, there are two approaches for comparing GGMs across any number of groups with either the posterior predictive distribution or with Bayesian hypothesis testing.

**BGGM** was primarily built for social-behvaioral scientists, where GGMs have become increasingly popular [@Epskamp2018]. For example, they have been used characterize the underlying structure of personality and psychopathology networks (i.e., non-zero partial correlations). These applications are distinct from other applications, such as gene co-expression [@Schafer2005] and functional connectivity applications  [@Das2017], in that high-dimensional data ($n < p$) is uncommon. This fact allows for extending inference beyond, say, $\ell_1$-regularized _point_ estimates of the partial correlations [@Friedman2008] and to seamlessly testing theorectical models with the Bayes factor. Additionally, by using Bayesian methodolgy, it is possible to explicitly gain evidence for 
conditional indepedence (i.e., the null hypothesis; $\mathcal{H}_0: \rho = 0$). 

Recently attention has shifted from estimating single networks to those from various sub-populations. **BGGM** also includes methods for comparing GGMs estimated from any number groups with either the predictive distribution or Bayesian hypothesis testing. Further, it is possible to test for invariances in graphical structures with the Bayes factor. This methodology is currently the only approach available  for comparing GGMs across any number of groups. Together, **BGGM** is a comprehensive R package that includes newly introduced methodology for Bayesian Gaussian graphical models.


## The Gaussian Graphical Model
The Gaussian graphical model captures conditional relationships \cite{lauritzen1996graphical} that are typically visualized to infer the underlying conditional (in)dependence structure \citep[i.e., the ``network";][]{Hojsgaard2012}. The undirected graph is $\mathcal{G} = (V, E)$, and includes a vertex set $V = \{1,...,p\}$ as well as an edge set $E \subset V \times V$. Let $\mathbf{y} = (y_1,...,y_p)^\top$ be a random vector indexed by the graphs vertices, of dimension $p$, that is assumed to follow a multivariate normal distribution $\mathcal{N}_p(\boldsymbol\mu, \boldsymbol\Sigma)$ and with a $p$ $\times$ $p$ positive definite covariance matrix $\boldsymbol\Sigma$. Without loss of information, the data is centered to have mean vector 0. Denote the precision matrix $\boldsymbol\Theta = \boldsymbol\Sigma^{-1}$. The graph is obtained from the off-diagonal elements $\theta_{ij} \in \boldsymbol\Theta_{ij}$. This is used to construct an adjacency matrix $A$ that follows

$$
\begin{equation}
A_{ij} = 
\begin{cases}
1, &  \text{if }\theta_{ij} \neq 0, \hspace{0.15 cm} 1 \leq i < j \leq p \\
0, & \text{otherwise},
\end{cases}
\end{equation}
$$
with $1 \leq i < j \leq p$ denoting the elements in the upper-triangular of the $p \times p$ matrix. Further, ($i,j$) $\in$ $E$ when the variables $i$ and $j$ are \emph{not} conditionally independent and  set to zero otherwise. Note that the edges are partial correlations ($\rho$) determined to be non-zero. These are computed directly from the precision matrix as

$$
\begin{align}
\label{eq:parcor}
    \rho_{ij} = \frac{- \theta_{ij}}{\sqrt{\theta_{ii}\theta_{jj}}},  \hspace{0.15 cm} 1 \leq i < j \leq p.
\end{align}
$$
Note that **BGGM** is organized around two general approaches for Bayesian inference--estimation X and hypothesis testing X. These partial correlations are explicitly used for the Bayes factor based approaches, whereas the  precision matrix is primarily used for the estimation based methods. 


### Estimation
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
Note that all functions are `S3` generics that include `summary` and `plot`. In this case, the generic is `estimate `. Further, by setting `analytic = F`, it is also possible to compare partial corrlations with `compare(.)` and nodewise predictability, including Bayesian $R^2$ X and (approximate) leave-one-out cross-validation X, can be computed with `predict(.)` (in- and out-of-sample). Bayesian $R^2$ can be computed with 

```r
fit_sample <- estimate(Y, analytic = F)

r2 <- predict(fit_analytic, 
              ci_width = .95, 
              measure = "R2") 
plot_1A <- plot(r2) +
            ggtitle("Network") +
            theme(title = element_text(size = 20))
```
The returned object is a ggplot X, which can then be futher customized.

### Hypothesis Testing
#### Exploratory
The following methods were introduced in Williams and Mulder (2019).  The Bayes factor is used for determing $E$--i.e.,

```r
# fit model
> fit_bf <- explore(Y, prior_sd = 0.35, 
+                   iter = 5000, cores = 2)
# select edge set
> E <- select(fit_bf, BF_cut = 3, 
+              alternative = "two.sided")
# summary
> summary(E, summarize = T, log = T, digits = 2)
# output

BGGM: Bayesian Gaussian Graphical Models 
--- 
Type: Hypothesis Testing 
Alternative: two.sided 
Bayes Factor: 3 
--- 
Call:
select.explore(x = fit_bf, BF_cut = 3, alternative = "two.sided")
--- 
Hypotheses: 
H0: rho = 0
H1: rho != 0 
--- 
Estimates: 
 
  edge post_mean post_sd BF 10
 1--2   -0.2405   0.018  76.4
 1--3   -0.1077   0.019  11.7
 2--3    0.2861   0.017 115.9
 1--4   -0.0075   0.019  -3.5
 2--4    0.1643   0.018  34.7
 3--4    0.1774   0.018  41.7
 1--5   -0.0086   0.019  -3.5
 2--5    0.1562   0.018  30.9
 3--5    0.3589   0.017 183.7
 4--5    0.1216   0.019  17.1
--- 
note: BF_10 is evidence in favor of H1
```
Note `prior_sd` governs the alternaitve hypothesis. It corresponds to the standard deviation of beta distribution (approximately) and it can be visualized with `hypothesis_plot(.)`. The graphical structure is visualized with `plot(E, type = "network")` (Figure 1 B). Further, one-sided hypothesis testing and an exhative approach ($\mathcal{H}_0: \rho_{ij} = 0$ vs. $\mathcal{H}_1: \rho_{ij} > 0$ vs. $\mathcal{H}_1: \rho_{ij} < 0$) is also an option (i.e., by changing `alternative`). The object `E` can then be plotted to visualize the network, for example, 

```r
plot_1B <- plot(E, type = "network", 
                labels = colnames(Y)) + 
                ggtitle("Network") +
                theme(title = element_text(size = 20))
```
Alternatively, a heatmap can be plotted by changing `type == "heatmap".`

#### Confirmatory

**BGGM** allows for comparing theoretical models, such as competing sets of order or equality constraints on multiple partial correlations in GGMs. This stands in contrast to data driven model selection that is commonplace in the GGM literature. These methods were introduced in X. For example, the following order contraint can be tested

$$
\begin{align}
\label{eq:order}
    \mathcal{H}_1&: (\rho_{1,2}, \rho_{1,3}) < 0 < \rho_{1,4} \\ \nonumber
    \mathcal{H}_c&: \text{``not} \hspace{.15 cm} \mathcal{H}_1 \text{"},
\end{align}
$$
which stats that $\rho_{1,2}$ and $\rho_{1,3}$ are predicted to be negative, whereas is $\rho_{1,4}$ is predicted to be positive. This is tested againts its compliment, $\text{``not} \hspace{.15 cm} \mathcal{H}_1 \text{"}$.

```{r}
# data
Y <- BGGM::bfi[,1:5]

# hypothesis
hypothesis <- c("(1--2, 1--3)  <  0 < (1--4)")

# fit model
test_order <-  confirm(Y, hypothesis = hypothesis, 
                       prior_sd = 0.5, 
                       iter = 50000, 
                       cores = 2)

# Output
BGGM: Bayesian Gaussian Graphical Models 
--- 
Type: Confirmatory Hypothesis Testing 
--- 
Hypotheses: 
                        
 H1 (1--2,1--3)<0<(1--4)
 Hc             'not H1'
--- 
Posterior prob: 
                 
 p(H1|Y) = 0.7906
 p(Hc|Y) = 0.2094
--- 
Bayes factor matrix: 
      H1    Hc
H1 1.000 0.265
Hc 3.775 1.000
--- 
note: equal hypothesis prior probabilities
```
Essentially any hypothesis can be specified, include inequality or equality contraints, as well as a combination of both. Further, several competeing hypothses can be tested.

### Comparing GGMs
