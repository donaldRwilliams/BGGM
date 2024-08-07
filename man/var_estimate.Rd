% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_estimate.R
\name{var_estimate}
\alias{var_estimate}
\title{VAR: Estimation}
\usage{
var_estimate(
  Y,
  rho_sd = sqrt(1/3),
  beta_sd = 1,
  iter = 5000,
  progress = TRUE,
  seed = NULL,
  ...
)
}
\arguments{
\item{Y}{Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).}

\item{rho_sd}{Numeric. Scale of the prior distribution for the partial correlations,
approximately the standard deviation of a beta distribution
(defaults to sqrt(1/3) as this results to delta = 2, and a uniform distribution across the partial correlations).}

\item{beta_sd}{Numeric. Standard deviation of the prior distribution for the regression coefficients
(defaults to 1). The prior is by default centered at zero and follows a normal distribution
\insertCite{@Equation 9, @sinay2014bayesian}{BGGM}}

\item{iter}{Number of iterations (posterior samples; defaults to 5000).}

\item{progress}{Logical. Should a progress bar be included (defaults to \code{TRUE}) ?}

\item{seed}{An integer for the random seed (defaults to 1).}

\item{...}{Currently ignored.}
}
\value{
An object of class \code{var_estimate} containing a lot of information that is
used for printing and plotting the results. For users of \strong{BGGM}, the following are the
useful objects:

\itemize{

\item \code{beta_mu} A matrix including the regression coefficients (posterior mean).

\item \code{pcor_mu} Partial correlation matrix (posterior mean).

\item \code{fit} A list including the posterior samples.

}
}
\description{
Estimate VAR(1) models by efficiently sampling from the posterior distribution. This
provides two graphical structures: (1) a network of undirected relations (the GGM, controlling for the
lagged predictors) and (2) a network of directed relations (the lagged coefficients). Note that
in the graphical modeling literature, this model is also known as a time series chain graphical model
\insertCite{abegaz2013sparse}{BGGM}.
}
\details{
Each time series in \code{Y} is standardized (mean  = 0; standard deviation = 1).
}
\note{
\strong{Regularization}:

A Bayesian ridge regression can be fitted by decreasing \code{beta_sd}
(e.g., \code{beta_sd = 0.25}). This could be advantageous for forecasting
(out-of-sample prediction) in particular.
}
\examples{
\donttest{
# data
Y <- subset(ifit, id == 1)[,-1]

# use alias (var_estimate also works)
fit <- var_estimate(Y, progress = FALSE)

fit

}
}
\references{
\insertAllCited{}
}
