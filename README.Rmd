# BGGM
Bayesian Gaussian Graphical Models

This package is described in Williams and Mulder (2019) and Williams (2018). The methods are seperated into two Bayesian approaches for inference: hypothesis testing and estimation. The former is described in Williams and Mulder (2018a), and allows for testing for the presence of edges with the Bayes factor. One-sided hypothesis testing is also possible. These methods can also provide evidence for the null hypothesis. There are extensions for confirmatory hypothesis testing in GGMs, that can include inequality or equality contraints on the partial correlation.

The estimation based method are described in Williams (2018). The methods offer advantages compared to classical methods, 
in that a measure of uncertainty is provided for all parameters. For example, each node has a distribution for
the variance explained. Measures of out-of-sample performance are also available. The model is selected with credible exclusion of zero.
<!---
#Comparing GGMs is also possible. These methods are described in Williams and Mulder (2019b). The comparison can be made with posterior #predictive loss functions for mulviarte normal distributions, or with the Bayes factor. Entire networks can be assess for similarities or #differences, or indidvudal nodes or edges can be tested. These methods can also be used to determine if network structures replicated #across two dataset, where evidence can be gained not only against but for hypothesized model.
-->

Download package here:

devtools::install_github("donaldRwilliams/BGGM")

Williams, D. R., & Mulder, J. (2019, January 14). Bayesian Hypothesis Testing for Gaussian Graphical Models: Conditional Independence and Order Constraints. https://doi.org/10.31234/osf.io/ypxd8

Williams, D. R. (2018, September 20). Bayesian Inference for Gaussian Graphical Models: Structure Learning, Explanation, and Prediction. https://doi.org/10.31234/osf.io/x8dpr
