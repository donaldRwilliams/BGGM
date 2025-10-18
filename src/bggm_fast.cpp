// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include <truncnorm.h>
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]

namespace {

const double kProbitCap = 8.0;

inline double Phi(double x) {
  return R::pnorm(x, 0.0, 1.0, TRUE, FALSE);
}

inline double PhiInv(double p) {
  return R::qnorm(p, 0.0, 1.0, TRUE, FALSE);
}

inline double StdNormPdf(double x) {
  return R::dnorm(x, 0.0, 1.0, FALSE);
}

inline double clamp_probit_bound(double value) {
  if (std::isfinite(value)) {
    if (value > kProbitCap) {
      return kProbitCap;
    }
    if (value < -kProbitCap) {
      return -kProbitCap;
    }
    return value;
  }
  return value > 0 ? kProbitCap : -kProbitCap;
}

double etruncnorm(double mu, double sigma, double lower, double upper) {
  double clipped_lower = clamp_probit_bound(lower);
  double clipped_upper = clamp_probit_bound(upper);
  double lower_bound = std::min(clipped_lower, clipped_upper);
  double upper_bound = std::max(clipped_lower, clipped_upper);
  double alpha = (clipped_lower - mu) / sigma;
  double beta = (clipped_upper - mu) / sigma;
  double Z = Phi(beta) - Phi(alpha);
  if (Z < 1e-12) {
    double midpoint = 0.5 * (clipped_lower + clipped_upper);
    if (!std::isfinite(midpoint)) {
      midpoint = mu;
    }
    return std::min(std::max(midpoint, lower_bound), upper_bound);
  }
  double numerator = StdNormPdf(alpha) - StdNormPdf(beta);
  double mean = mu + sigma * numerator / Z;
  mean = std::min(std::max(mean, lower_bound), upper_bound);
  return mean;
}

double truncated_cdf_diff(double mu, double sigma, double lower, double upper) {
  double clipped_lower = clamp_probit_bound(lower);
  double clipped_upper = clamp_probit_bound(upper);
  double lower_cdf = Phi((clipped_lower - mu) / sigma);
  double upper_cdf = Phi((clipped_upper - mu) / sigma);
  double diff = upper_cdf - lower_cdf;
  return std::max(diff, 1e-12);
}

double rtnorm(double mu, double sigma, double lower, double upper) {
  double clipped_lower = clamp_probit_bound(lower);
  double clipped_upper = clamp_probit_bound(upper);
  double lower_cdf = Phi((clipped_lower - mu) / sigma);
  double upper_cdf = Phi((clipped_upper - mu) / sigma);
  double diff = upper_cdf - lower_cdf;
  if (diff < 1e-12) {
    return etruncnorm(mu, sigma, clipped_lower, clipped_upper);
  }
  double u = R::runif(lower_cdf, upper_cdf);
  return R::qnorm(u, mu, sigma, TRUE, FALSE);
}

arma::rowvec init_thresholds_from_data(const arma::vec& y, std::size_t num_levels) {
  arma::rowvec thresholds(num_levels + 1, arma::fill::zeros);
  thresholds[0] = -arma::datum::inf;
  thresholds[num_levels] = arma::datum::inf;
  double n = static_cast<double>(y.n_elem);
  double min_step = 1e-6;
  for (std::size_t c = 1; c < num_levels; ++c) {
    double cum_prob = arma::accu(y <= static_cast<double>(c)) / n;
    cum_prob = std::min(std::max(cum_prob, min_step), 1.0 - min_step);
    double proposed = PhiInv(cum_prob);
    if (c > 1 && proposed <= thresholds[c - 1]) {
      thresholds[c] = thresholds[c - 1] + min_step;
    } else {
      thresholds[c] = proposed;
    }
  }
  return thresholds;
}

void recenter_column_and_thresholds(arma::cube& Z,
                                    int col_index,
                                    arma::mat& thresholds,
                                    arma::mat* candidate = NULL) {
  arma::vec column = Z.slice(0).col(col_index);
  double mean_val = arma::mean(column);
  Z.slice(0).col(col_index) -= mean_val;
  for (arma::uword idx = 0; idx < thresholds.n_cols; ++idx) {
    double& value = thresholds(col_index, idx);
    if (std::isfinite(value)) {
      value -= mean_val;
    }
  }
  if (candidate != NULL) {
    for (arma::uword idx = 0; idx < candidate->n_cols; ++idx) {
      double& value = (*candidate)(col_index, idx);
      if (std::isfinite(value)) {
        value -= mean_val;
      }
    }
  }
}

void init_Z_from_trunc_means(const arma::mat& Y,
                             arma::cube& Z,
                             arma::mat& thresholds) {
  int n = Y.n_rows;
  int k = Y.n_cols;
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row < n; ++row) {
      int category = static_cast<int>(Y(row, col));
      category = std::max(category, 1);
      category = std::min(category, static_cast<int>(thresholds.n_cols) - 1);
      double lower = thresholds(col, category - 1);
      double upper = thresholds(col, category);
      Z(row, col, 0) = etruncnorm(0.0, 1.0, lower, upper);
    }
    recenter_column_and_thresholds(Z, col, thresholds);
  }
}

}  // namespace

// mean of 3d array
// [[Rcpp::export]]
arma::mat mean_array(arma::cube x){
  return mean(x, 2);
}

// R quantile type = 1
// [[Rcpp::export]]
double quantile_type_1(arma::vec x, double prob){

  arma::mat sort_x = sort(x.elem(find_finite(x)));

  int n = sort_x.n_rows;

  float nppm =  n * prob;

  float j = floor(nppm);

  //float h = 0;

  float qs = 0;

  // if(nppm > j){
  //
  //   float h = 1;
  //
  // } else {
  //
  //   float h = 0;
  //
  // }

  arma::mat x_1(2, 1);
  arma::mat x_n(2, 1);

  x_1.col(0).row(0) = sort_x(0);
  x_1.col(0).row(1) = sort_x(0);

  x_n.col(0).row(0) = sort_x(n - 1);
  x_n.col(0).row(1) = sort_x(n - 1);

  arma::mat join_x = join_vert(x_1, sort_x, x_n);

  if(nppm > j){

    qs = join_x(j + 3);

  } else {

    qs = join_x(j + 2);

  }

  return qs;

}


// [[Rcpp::export]]
arma::mat Sigma_i_not_i(arma::mat x, int index) {
  arma::mat sub_x = x.row(index);
  sub_x.shed_col(index);
  return(sub_x);
}

// [[Rcpp::export]]
arma::vec select_col(arma::mat x, int index){
  arma::vec z = x.col(index);
  return(z);
}

// [[Rcpp::export]]
arma::mat select_row(arma::mat x, int index){
  arma::mat z = x.row(index);
  return(z);
}

// [[Rcpp::export]]
arma::mat remove_row(arma::mat x, int which){
  x.shed_row(which);
  return(x);
}

// [[Rcpp::export]]
arma::mat remove_col(arma::mat x, int index){
  x.shed_col(index);
  return(x);
}


// Hoff, P. D. (2009). A first course in Bayesian statistical
// methods (Vol. 580). New York: Springer.
// pp 105-123

// note: `internal` is a simplified version of missing_gaussian
// that seemed to be faster when used within the Gibbs sampler

// [[Rcpp::export]]
Rcpp::List internal_missing_gaussian(arma::mat Y,
                            arma::mat Y_missing,
                            arma::mat Sigma,
                            int iter_missing) {
  int p = Y.n_cols;
  int n = Y.n_rows;

  arma::uvec index = find(Y_missing == 1);

  int n_na = index.n_elem;

  arma::mat ppc_missing(iter_missing, n_na, arma::fill::zeros);

  for(int s = 0; s < iter_missing; ++s){

    for(int j = 0; j < p; ++j){

      arma::vec Y_j = Y_missing.col(j);

      double check_na = sum(Y_j);

      if(check_na == 0){
        continue;
      }

      arma::uvec  index_j = find(Y_missing.col(j) == 1);
      int  n_missing = index_j.n_elem;

      arma::mat beta_j = Sigma_i_not_i(Sigma, j) * inv(remove_row(remove_col(Sigma, j), j));
      arma::mat  sd_j = sqrt(select_row(Sigma, j).col(j) - Sigma_i_not_i(Sigma, j) *
      inv(remove_row(remove_col(Sigma, j), j)) * Sigma_i_not_i(Sigma, j).t());
      arma::vec pred = remove_col(Y,j) * beta_j.t();
      arma::vec pred_miss = pred(index_j);

      for(int i = 0; i < n_missing; ++i){
        arma::vec ppc_i = Rcpp::rnorm(1,  pred(index_j[i]), arma::as_scalar(sd_j));
        Y.col(j).row(index_j[i]) = arma::as_scalar(ppc_i);
      }
    }

    arma::mat S_Y = Y.t() * Y;
    arma::mat Theta = wishrnd(inv(S_Y),   (n - 1));
    Sigma = inv(Theta);
    ppc_missing.row(s) = Y.elem(index).t();
  }
  Rcpp::List ret;
  ret["Y"] = Y;
  ret["ppc_missing"] = ppc_missing;
  return  ret;
}



// Hoff, P. D. (2009). A first course in Bayesian statistical
// methods (Vol. 580). New York: Springer.
// pp 105-123

// [[Rcpp::export]]
Rcpp::List missing_gaussian(arma::mat Y,
                            arma::mat Y_missing,
                            arma::mat Sigma,
                            int iter_missing,
                            bool progress_impute,
                            bool store_all,
                            float lambda) {


  // progress
  Progress  pr(iter_missing, progress_impute);

  int p = Y.n_cols;
  int n = Y.n_rows;

  arma::uvec index = find(Y_missing == 1);

  int n_na = index.n_elem;

  arma::mat I_p(p, p, arma::fill::eye);

  // store posterior predictive distribution for missing values
  arma::mat ppd_missing(iter_missing, n_na, arma::fill::zeros);

  // store all imputed data sets
  arma::cube Y_all(n, p, iter_missing, arma::fill::zeros);

  for(int s = 0; s < iter_missing; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    for(int j = 0; j < p; ++j){

      arma::vec Y_j = Y_missing.col(j);

      double check_na = sum(Y_j);

      if(check_na == 0){
        continue;
      }

      arma::uvec  index_j = find(Y_missing.col(j) == 1);

      int n_missing = index_j.n_elem;

      arma::mat beta_j = Sigma_i_not_i(Sigma, j) * inv(remove_row(remove_col(Sigma, j), j));

      arma::mat  sd_j = sqrt(select_row(Sigma, j).col(j) - Sigma_i_not_i(Sigma, j) *
        inv(remove_row(remove_col(Sigma, j), j)) * Sigma_i_not_i(Sigma, j).t());

      arma::vec pred = remove_col(Y,j) * beta_j.t();

      arma::vec pred_miss = pred(index_j);

      for(int i = 0; i < n_missing; ++i){

        arma::vec ppd_i = Rcpp::rnorm(1,  pred(index_j[i]),
                                      arma::as_scalar(sd_j));

        Y.col(j).row(index_j[i]) = arma::as_scalar(ppd_i);

      }

    }

    arma::mat S_Y = Y.t() * Y;
    arma::mat Theta = wishrnd(inv(S_Y + I_p * lambda), n + lambda);
    Sigma = inv(Theta);

    if(store_all){
      Y_all.slice(s) = Y;
    }
  }

  Rcpp::List ret;
  ret["Y_all"] = Y_all;
  return ret;
}




// matrix F continous sampler
// Williams, D. R., & Mulder, J. (2019). Bayesian hypothesis testing for Gaussian
// graphical models:  Conditional independence and order constraints.

// [[Rcpp::export]]
Rcpp::List Theta_continuous(arma::mat Y,
                            int iter,
                            float delta,
                            float epsilon,
                            int prior_only,
                            int explore,
                            arma::mat start,
                            bool progress,
                            bool impute,
                            arma::mat Y_missing) {



  // note p changed to k to be consistent
  //with the multivariate regression samplers
  Progress  p(iter, progress);

  // number of rows
  float n = Y.n_rows;

  int k = 1;

  // sample prior
  if(prior_only == 1){

    if(explore == 1){

      k = 3;

      } else {

      k = Y.n_cols;
        }

      } else {

    // number of columns
    k = Y.n_cols;

        }

  // k by k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  int nu = 1/ epsilon;

  // // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;

  // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  arma::mat B(epsilon * I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  arma::cube Sigma(k, k, 1, arma::fill::zeros);

  // starting value
  Psi.slice(0).fill(arma::fill::eye);

  Theta.slice(0) = start;

  arma::mat S_Y(Y.t() * Y);

  arma::uvec index = find(Y_missing == 1);

  int n_na = index.n_elem;

  arma::mat ppd_missing(iter, n_na, arma::fill::zeros);

  float iter_missing = 1;

  for(int  s = 0; s < iter; ++s){

    p.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    if(prior_only == 1){

      Psi.slice(0) = wishrnd(I_k * epsilon, nu);

      // sample Theta
      Sigma.slice(0) =   wishrnd(inv(Psi.slice(0)),   k - 1 + delta);

      // Sigma
      Theta.slice(0) = inv(Sigma.slice(0));


    } else {

      Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

      // sample Theta
      Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

    }

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    // store posterior samples
    pcors_mcmc.slice(s) =  -(pcors - I_k);


    if(impute){

      Rcpp::List ppd_impute = internal_missing_gaussian(Y, Y_missing,
                                               inv(Theta.slice(0)),
                                               iter_missing);
      // imputed Y
      arma::mat Y = ppd_impute["Y"];

      // scatter matrix
      S_Y = Y.t() * Y;

      // store missing values
      ppd_missing.row(s) = Y.elem(index).t();
    }

  }

  arma::cube fisher_z = atanh(pcors_mcmc);

  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter - 50), 2);

  arma::mat  ppd_mean = mean(ppd_missing, 0).t();

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] =  pcor_mat;
  ret["fisher_z"] = fisher_z;
  ret["ppd_mean"] = ppd_mean;
  return ret;
}


// [[Rcpp::export]]
Rcpp::List sample_prior(arma::mat Y,
                            int iter,
                            float delta,
                            float epsilon,
                            int prior_only,
                            int explore,
                            bool progress) {

  // note p changed to k to be consistent
  //with the multivariate regression samplers

  Progress  p(iter, progress);

  // number of rows
  float n = Y.n_rows;

  int k = 1;

  // sample prior
  if(prior_only == 1){

    if(explore == 1){

      k = 3;

    } else {

      k = Y.n_cols;

    }

  } else {

    // number of columns
    k = Y.n_cols;

  }

  // k by k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  int nu = 1 / epsilon;
  // // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;
  //
  // // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  arma::mat B(epsilon * I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);

  // covariance matrix
  arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);
  arma::cube Sigma(k, k, 1, arma::fill::zeros);

  // starting value
  Sigma.slice(0).fill(arma::fill::eye);
  Psi.slice(0).fill(arma::fill::eye);
  Theta.slice(0).fill(arma::fill::eye);

  arma::mat S_Y(Y.t() * Y);

  for(int  s = 0; s < iter; ++s){

    p.increment();

    if(s % 250 == 0){

      Rcpp::checkUserInterrupt();

    }


    if(prior_only == 1){

      Psi.slice(0) = wishrnd(I_k * epsilon, nu);

      // sample Theta
      Sigma.slice(0) =   wishrnd(inv(Psi.slice(0)),   k - 1 + delta);

      // Sigma
      Theta.slice(0) = inv(Sigma.slice(0));


    } else {

      Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

      // sample Theta
      Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

      // Sigma
      Sigma.slice(0) = inv(Theta.slice(0));

    }

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    // store posterior samples
    pcors_mcmc.slice(s) =  -(pcors - I_k);

    }

  arma::cube fisher_z = atanh(pcors_mcmc);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["fisher_z"] = fisher_z;
  return ret;
}


// [[Rcpp::export]]
Rcpp::List mv_continuous(arma::mat Y,
                          arma::mat X,
                          float delta,
                          float epsilon,
                          int iter,
                          arma::mat start,
                          bool progress){


  // progress
  Progress  pr(iter, progress);

  // number of rows
  int n = Y.n_rows;

  // number of dependent variables
  int k = Y.n_cols;

  // number of predictors
  int p = X.n_cols;

  int nu = 1/ epsilon;

  // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1;

  // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // k * k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  // p * p identity mat
  arma::mat  I_p(p, p, arma::fill::eye);

  // scatter matrix X' * X
  arma::mat S_X(X.t() * X + I_p * 0.000001);

  // inv S_X
  arma::mat Sinv_X(inv(S_X));

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  // scatter matrix dependent variables
  arma::mat S_Y(k, k, arma::fill::zeros);
  arma::mat B(epsilon*I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);

  // covariance matrix
  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  // coefficients
  arma::mat beta(p, k, arma::fill::zeros);
  arma::cube beta_mcmc(p, k, iter,  arma::fill::zeros);

  // starting value
  Sigma.slice(0) = inv(start);
  Psi.slice(0).fill(arma::fill::eye);
  Theta.slice(0) = start;

  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    // draw coefficients
    beta = reshape(mvnrnd(reshape(Sinv_X * X.t() * Y, k*p , 1),
                          kron(Sigma.slice(0), Sinv_X)),
                          p, k);

    // scatter matrix
    S_Y = Y.t() * Y + I_k - beta.t() * S_X * beta;

    // sample Psi
    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

    // sample Theta
    Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

    // Sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));


    beta_mcmc.slice(s) = beta;
    pcors_mcmc.slice(s) =  -(pcors - I_k);
  }

  arma::cube fisher_z = atanh(pcors_mcmc);

  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter - 50), 2);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] =  pcor_mat;
  ret["beta"] = beta_mcmc;
  ret["fisher_z"] = fisher_z;
  return ret;
}

// [[Rcpp::export]]
Rcpp::List trunc_mvn(arma::mat mu,
                     arma::mat rinv,
                     arma::mat z,
                     arma::mat y,
                     arma::rowvec cutpoints){

  // adapted from Aline Talhouk (matlab)
  // mu: beta hat (p x 1)
  // rinv: matrix-F precision matrix
  // z: latent data (p x 1)
  // y: observed binary data (for now)
  // cutpoint: thresholds for truncati

  // number of columns
  int T = rinv.n_cols;

  // number of cutpoints
  int n_cuts = cutpoints.n_cols;

  // substrat mean
  arma::mat zt(z - mu);

  // cutpoints (a: low, b: upper)
  arma::mat a(1, n_cuts - 1); a.fill(0);
  arma::mat b(1, n_cuts - 1); a.fill(0);

  for(int i = 0; i < (n_cuts - 1); ++i){

    a.col(i) = cutpoints(i);

    b.col(i) = cutpoints(i + 1);

  }

  // match cutpoints to T
  arma::mat a2(T, 1);
  arma::mat b2(T, 1);

  for(int i = 0; i < T; ++i){

    a2.row(i) = a.col(arma::as_scalar(y.row(i)));

    b2.row(i) = b.col(arma::as_scalar(y.row(i)));

  }

  // alpha mat
  arma::mat alpha(a2 - mu);
  // beta mat
  arma::mat beta(b2 - mu);

  // h2
  arma::mat h2(1 / rinv.diag());

  // temporary holders to make c
  arma::mat diag1(T, T);
  arma::mat diag2 = rinv.diag();

  for(int i = 0; i < T; ++i){

    diag1.row(i) = arma::repmat(diag2.row(i), 1, T);

  }

  // c
  arma::mat c(- rinv / diag1);
  for(int i = 0; i < T; ++i){
    c.row(i).col(i) = 0;
  }

  // term
  arma::mat  term(T, 1, arma::fill::zeros);

  // upper bound
  arma::mat  lb(T, 1, arma::fill::zeros);

  // lower bound
  arma::mat  ub(T, 1, arma::fill::zeros);

  // epsilon
  // arma:vec  eps(1, fill::zeros);

  for(int i = 0; i < T; ++i){

    arma::mat term(c * zt);

    arma::mat lb((alpha - repmat(term.row(i), T, 1)) / sqrt(h2(i)));

    arma::mat ub((beta - repmat(term.row(i), T, 1)) / sqrt(h2(i)));

    arma::vec  eps(rtruncnorm(1, 0,  1,  lb(i),  ub(i)));

    zt.row(i) = term.row(i)  + sqrt(h2(i)) *  eps;

  }

  arma::mat zr(zt + mu);

  // return
  Rcpp::List ret;
  ret["zr"] = zr;
  return  ret;
}

// binary sampler
// [[Rcpp::export]]
Rcpp::List mv_binary(arma::mat Y,
		     arma::mat X,
		     float delta,
		     float epsilon,
		     int iter,
		     float beta_prior,
		     arma::rowvec cutpoints,
		     arma::mat start,
		     bool progress){

  // Y: data matrix (n * k)
  // X: predictors (n * p) (blank for "network")
  // iter: number of samples
  // cutpoints: truncation points
  // delta: hyperparameter

  // progress
  Progress  pr(iter, progress);

  // dependent variables
  int k = Y.n_cols;

  // sample size
  int n = Y.n_rows;

  // number of predictors
  int p = X.n_cols;

  // arma::rowvec ct = cutpoints;

  // int epsilon1 = epsilon;

  int nu = 1 / epsilon;
  // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;

  // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // k * k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  // p * p identity mat
  arma::mat  I_p(p, p, arma::fill::eye);

  // scatter matrix X' * X
  arma::mat S_X(X.t() * X + I_p * beta_prior);

  // inv S_X
  arma::mat Sinv_X(inv(S_X));

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  // scatter matrix dependent variables
  arma::mat S_Y(k, k, arma::fill::zeros);
  arma::mat B(epsilon * I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  // arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);

  // covariance matrix
  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  // arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  // latent data
  arma::cube z0(n, k, 1,  arma::fill::zeros);

  // expanded latent data
  arma::mat w(n, k, arma::fill::zeros);

  // conditonal data
  arma::cube Xbhat(n, k,  1, arma::fill::zeros);

  // Rinv update
  arma::cube Rinv(k, k, 1, arma::fill::zeros);

  arma::cube R(k, k, 1, arma::fill::zeros);

  // coefficients
  arma::mat M(p, k, arma::fill::zeros);

  // expanded coefs
  arma::mat beta(p, k, arma::fill::zeros);
  arma::cube beta_mcmc(p, k, iter,  arma::fill::zeros);

  // draw coefs conditional on w
  arma::mat gamma(p, k, arma::fill::zeros);

  arma::cube Dinv(k, k, 1, arma::fill::zeros);
  arma::mat D(k, k,  arma::fill::eye);

  // starting values
  Sigma.slice(0) = inv(start);
  Theta.slice(0) = start;
  Psi.slice(0).fill(arma::fill::eye);
  Dinv.slice(0).fill(arma::fill::eye);
  Rinv.slice(0).fill(arma::fill::eye);
  R.slice(0).fill(arma::fill::eye);

  arma::mat ss(1, 1);
  arma::mat mm(n,1);

  // start sampling
  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }


    for(int i = 0; i < k; ++i){

      mm = Xbhat.slice(0).col(i).t() +
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        (remove_col(z0.slice(0), i).t() - remove_col(Xbhat.slice(0), i).t());

      ss = select_row(R.slice(0), i).col(i) -
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        Sigma_i_not_i(R.slice(0), i).t();


      for(int j = 0 ; j < n; ++j){

        arma::vec temp = Y.col(i).row(j);

        if(temp(0) == 0){


          arma::vec temp_j = rtruncnorm(1,   mm(j),  sqrt(ss(0)),  -arma::datum::inf,  0);

          z0.slice(0).col(i).row(j) =   temp_j;


        } else {



          arma::vec temp_j = rtruncnorm(1,   mm(j),  sqrt(ss(0)),  0, arma::datum::inf);

          z0.slice(0).col(i).row(j) =   temp_j;

        }
      }
    }

    // D matrix
    for(int i = 0; i < k; ++i){

      D.row(i).col(i) = sqrt(1 / R::rgamma((delta + k - 1) / 2,
                   2 / arma::as_scalar(Rinv.slice(0).row(i).col(i))));
    }

    w = z0.slice(0) * D;
    
    M  = Sinv_X * X.t() * w;

    gamma = reshape(mvnrnd(reshape(M, k * p , 1),
                           kron(Sigma.slice(0), Sinv_X)),
                           p, k);

    beta  = gamma * Dinv.slice(0);

    Xbhat.slice(0) = X * beta;

    S_Y =   w.t() * w + I_k - M.t() * S_X * M;

    // sample Psi
    // Debugging:
    // Declare the bmpinv variable
    // arma::mat bmpinv;
    // bmpinv = inv(BMPinv + Theta.slice(0));
    // if (!bmpinv.is_sympd()) {
    //   // Print the bmpinv matrix for debugging
    //   bmpinv.print("bmpinv matrix:");
    //   Rcpp::stop("bmpinv matrix is not symmetric positive definite.");
    // }
    // END Debug

    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

    // sample Theta
    // Debugging:
    // Declare the psisl variable
    // arma::mat psisl;
    // psisl = inv(Psi.slice(0) + S_Y);
    // if (!psisl.is_sympd()) {
    //   // Print the matrix for debugging
    //   Psi.slice(0).print("Psi.slice(0):" );
    //   S_Y.print("S_Y:");
    //   w.print("w:" );
    //   I_k.print("I_k:");
    //   M.print("M:");
    //   S_X.print("S_X:");
    //   psisl.print("psisl matrix:");
    //   Rcpp::stop("psisl matrix is not symmetric positive definite.");
    // }
    // END Debug

    Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

    // Sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));


    Dinv.slice(0)  = inv(diagmat(sqrt(Sigma.slice(0).diag())));

    Rinv.slice(0)   = inv(cors);

    R.slice(0) = cors;

    beta_mcmc.slice(s) =reshape(beta, p,k);
    pcors_mcmc.slice(s) =  -(pcors - I_k);

  }

  arma::cube fisher_z = atanh(pcors_mcmc);
  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter - 50), 2);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] = pcor_mat;
  ret["beta"] = beta_mcmc;
  ret["fisher_z"] = fisher_z;
  return  ret;
}


// ordinal sampler
// [[Rcpp::export]]
Rcpp::List mv_ordinal_cowles(arma::mat Y,
                           arma::mat X,
                           float delta,
                           float epsilon,
                           int iter, float MH) {

  // number of rows
  float n = Y.n_rows;

  // number of columns
  int k = Y.n_cols;

  // number of predictors
  int p = X.n_cols;

  int nu = 1/ epsilon;
  // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;

  // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // k by k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  // p by p identity mat
  arma::mat  I_p(p, p, arma::fill::eye);

  // scatter matrix X' * X
  arma::mat S_X(X.t() * X + I_p * 0.000001);

  // inv S_X
  arma::mat Sinv_X(inv(S_X));

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  // scatter matrix dependent variables
  arma::mat S_Y(k, k, arma::fill::zeros);
  arma::mat B(epsilon * I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);


  // covariance matrix
  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  // coefficients
  arma::mat beta(p, k, arma::fill::zeros);
  arma::cube beta_mcmc(p, k, iter,  arma::fill::zeros);

  // latent update
  arma::cube z0(n, k, 1,  arma::fill::zeros);

  // expanded latent data
  arma::mat w(n, k, arma::fill::zeros);

  // x hat
  arma::cube Xbhat(n, k,  1, arma::fill::zeros);

  // Rinv update
  arma::cube Rinv(k, k, 1, arma::fill::zeros);

  // coefficients
  arma::mat M(p, k, arma::fill::zeros);

  // Rinv update
  arma::cube R(k, k, 1, arma::fill::zeros);


  // Dinv update
  arma::cube Dinv(k, k, 1, arma::fill::zeros);

  // Dinv update
  arma::mat D(k, k, arma::fill::eye);

  // ordinal levels
  arma::vec n_levels = unique(Y.col(0));

  // MH step acceptence
  arma::vec acc(k, arma::fill::zeros);

  Rinv.slice(0).fill(arma::fill::eye);
  R.slice(0).fill(arma::fill::eye);
  Dinv.slice(0).fill(arma::fill::eye);
  Psi.slice(0).fill(arma::fill::eye);
  Sigma.slice(0).fill(arma::fill::eye);

  // range
  Rcpp::Range thresh_sampled =  Rcpp::seq(min(Y.col(0)) + 1, max(Y.col(0)) - 1);

  // threshold holders
  arma::cube  thresh_mat(k, n_levels.size() + 1, 1, arma::fill::zeros);
  arma::cube c_thresh_mat(k, n_levels.size() + 1, 1, arma::fill::zeros);

  // set ends to -Inf and Inf
  thresh_mat.slice(0).col(0).fill(-arma::datum::inf);
  thresh_mat.slice(0).col(n_levels.size()).fill(arma::datum::inf);
  c_thresh_mat.slice(0).col(0).fill(-arma::datum::inf);
  c_thresh_mat.slice(0).col(n_levels.size()).fill(arma::datum::inf);

  // store thresholds
  arma::cube thresh_mcmc(k, n_levels.size() + 1, iter);

  // MH R ratio
  arma::mat R_numer(n, k);
  arma::mat R_denom(n, k);

  arma::mat mm(n,1);
  arma::mat ss(1,1);


  // draw coefs conditional on w
  arma::mat gamma(p, k, arma::fill::zeros);

  arma::mat current_thresh = thresh_mat.slice(0);
  arma::mat candidate_thresh = c_thresh_mat.slice(0);

  // starting thresholds
  for (int j = 0; j < k; ++j) {
    arma::rowvec tau = init_thresholds_from_data(Y.col(j), n_levels.size());
    current_thresh.row(j) = tau;
  }
  candidate_thresh = current_thresh;

  init_Z_from_trunc_means(Y, z0, current_thresh);
  candidate_thresh = current_thresh;
  thresh_mat.slice(0) = current_thresh;
  c_thresh_mat.slice(0) = candidate_thresh;
  thresh_mcmc.slice(0) = current_thresh;

  // start sampling
  for(int s = 0; s < iter; ++s){

    Rcpp::checkUserInterrupt();

    candidate_thresh = current_thresh;

    // multivariate data
    for(int i = 0; i < k; ++i){

      arma::mat mm = Xbhat.slice(0).col(i).t() +
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        (remove_col(z0.slice(0), i).t() - remove_col(Xbhat.slice(0), i).t());

      arma::vec ss = select_row(R.slice(0), i).col(i) -
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        Sigma_i_not_i(R.slice(0), i).t();

      // generate latent data
      double sd = std::sqrt(ss(0));
      for(int j = 0; j < n; ++j){
        int category = static_cast<int>(Y.col(i)[j]);
        double lower = current_thresh(i , category - 1);
        double upper = current_thresh(i , category);
        z0.slice(0).row(j).col(i) = rtnorm(mm(j), sd, lower, upper);
      }

      recenter_column_and_thresholds(z0, i, current_thresh);
      candidate_thresh.row(i) = current_thresh.row(i);
    }

    for(int i = 0; i < max(Y.col(0)) - 1; ++i){
      for(int j = 0; j < k ; ++j){
        double location = current_thresh(j , thresh_sampled[i]);
        candidate_thresh(j, i+2) = rtnorm(location, MH, 0.0, arma::datum::inf);
      }
    }

    for(int j = 0; j < n; ++j){
      for(int i = 0; i < k; ++i){
        int category = static_cast<int>(Y.col(i)[j]);
        double mu = mm(j);
        double sd = std::sqrt(ss(0));
        R_numer(j, i) = truncated_cdf_diff(mu, sd,
                                           candidate_thresh(i , category - 1),
                                           candidate_thresh(i , category));

        R_denom(j, i) = truncated_cdf_diff(mu, sd,
                                           current_thresh(i , category - 1),
                                           current_thresh(i , category));
      }
    }

    // remove nan
    R_numer.replace(arma::datum::nan, 0);
    R_denom.replace(arma::datum::nan, 0);

    // block update for each variable
    for(int j = 0; j < k; j ++){
      if(accu(log(R_numer.col(j) / R_denom.col(j))) > log(R::runif(0, 1))){
        // acceptence
        acc(j) = acc(j) + 1;
        // block update
        current_thresh.row(j) = candidate_thresh.row(j);
      }
    }

    for(int i = 0; i < k; ++i){
      D.row(i).col(i) = sqrt(1 / R::rgamma((deltaMP + k - 1) / 2,
                   2 / arma::as_scalar(Rinv.slice(0).row(i).col(i))));
    }

    // expand latent data
    w = z0.slice(0) * D;

    // coefficients
    M = Sinv_X * X.t() * w;

    // gamma
    gamma = reshape(mvnrnd(reshape(M, k * p , 1),
                           kron(Sigma.slice(0), Sinv_X)),
                           p, k);

    // expanded coefficients
    beta  = gamma * Dinv.slice(0);

    // update Yhat
    Xbhat.slice(0) = X * beta;

    // error scatter matrix
    S_Y =   w.t() * w + I_k - M.t() * S_X * M;

    // } while (det(S_Y) < 0);

    // sample Psi
    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

    // sample Theta
    Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

    // sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    // update Dinv
    Dinv.slice(0)  = inv(diagmat(sqrt(Sigma.slice(0).diag())));

    // inverse correlation matrix
    Rinv.slice(0)   = inv(cors);

    // update correlation matrix
    R.slice(0) = cors;

    beta_mcmc.slice(s) =reshape(beta, p,k);
    pcors_mcmc.slice(s) =  -(pcors - I_k);
    cors_mcmc.slice(s) =  cors;
    Sigma_mcmc.slice(s) = Sigma.slice(0);
    Theta_mcmc.slice(s) = Theta.slice(0);
    thresh_mat.slice(0) = current_thresh;
    c_thresh_mat.slice(0) = candidate_thresh;
    thresh_mcmc.slice(s) = current_thresh;

  }

  ////////////////////////////////////////////

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["cors"] =  cors_mcmc;
  ret["beta"] = beta_mcmc;
  ret["Theta"] = Theta_mcmc;
  ret["Sigma"] = Sigma_mcmc;
  ret["acc"]  = acc;
  ret["thresh"]  = thresh_mcmc;
  return  ret;
}

// ordinal sampler customary
// [[Rcpp::export]]
Rcpp::List mv_ordinal_albert(arma::mat Y,
                          arma::mat X,
                          int iter,
                          float delta,
                          float epsilon,
                          int K,
                          arma::mat start,
                          bool progress
                          ){


  Progress  pr(iter, progress);

  // number of rows
  float n = Y.n_rows;

  // number of columns
  int k = Y.n_cols;

  // number of predictors
  int p = X.n_cols;


  // ordinal levels
  // int K = unique(Y.col(0));

  int nu = 1/ epsilon;
  // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;

  // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // k by k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  // p by p identity mat
  arma::mat  I_p(p, p, arma::fill::eye);

  // scatter matrix X' * X
  arma::mat S_X(X.t() * X + I_p * 0.001);

  // inv S_X
  arma::mat Sinv_X(inv(S_X));

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  // scatter matrix dependent variables
  arma::mat S_Y(k, k, arma::fill::zeros);
  arma::mat B(epsilon * I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  // arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);


  // covariance matrix
  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  // arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  // coefficients
  arma::mat beta(p, k, arma::fill::zeros);
  arma::cube beta_mcmc(p, k, iter,  arma::fill::zeros);

  // latent update
  arma::cube z0(n, k, 1,  arma::fill::zeros);

  // expanded latent data
  arma::mat w(n, k, arma::fill::zeros);

  // x hat
  arma::cube Xbhat(n, k,  1, arma::fill::zeros);

  // Rinv update
  arma::cube Rinv(k, k, 1, arma::fill::zeros);

  // coefficients
  arma::mat M(p, k, arma::fill::zeros);

  // Rinv update
  arma::cube R(k, k, 1, arma::fill::zeros);


  // Dinv update
  arma::cube Dinv(k, k, 1, arma::fill::zeros);

  // Dinv update
  arma::mat D(k, k, arma::fill::eye);

  // ordinal levels
  arma::vec n_levels = unique(Y.col(0));

  Rinv.slice(0).fill(arma::fill::eye);
  R.slice(0).fill(arma::fill::eye);
  Dinv.slice(0).fill(arma::fill::eye);
  Psi.slice(0).fill(arma::fill::eye);
  Sigma.slice(0) = inv(start);
  Theta.slice(0) = start;
  // draw coefs conditional on w
  arma::mat gamma(p, k, arma::fill::zeros);

  arma::cube thresh(iter, K+1, k, arma::fill::zeros);
  arma::mat current_thresh(k, K+1, arma::fill::zeros);

  for (int j = 0; j < k; ++j) {
    current_thresh.row(j) = init_thresholds_from_data(Y.col(j), K);
    current_thresh(j, 0) = -arma::datum::inf;
    current_thresh(j, K) = arma::datum::inf;
    thresh.slice(j).col(0).fill(-arma::datum::inf);
    thresh.slice(j).col(K).fill(arma::datum::inf);
    thresh.slice(j).row(0) = current_thresh.row(j);
  }

  init_Z_from_trunc_means(Y, z0, current_thresh);

  for (int j = 0; j < k; ++j) {
    thresh.slice(j).row(0) = current_thresh.row(j);
  }

  arma::mat mm(n,1);
  arma::mat ss(1,1);

  for(int s = 1; s < iter; ++s ){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    for(int i = 0; i < k; ++i){

      mm = Xbhat.slice(0).col(i).t() +
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        (remove_col(z0.slice(0), i).t() - remove_col(Xbhat.slice(0), i).t());

      ss = select_row(R.slice(0), i).col(i) -
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        Sigma_i_not_i(R.slice(0), i).t();

      double sd = std::sqrt(ss(0));

      if (s > 1) {
        for (int level = 2; level < K; ++level) {
          double lower_candidate = current_thresh(i, level - 1);
          double upper_candidate = current_thresh(i, level + 1);
          arma::uvec lower_idx = arma::find(Y.col(i) == level);
          if (!lower_idx.is_empty()) {
            arma::vec z_col = z0.slice(0).col(i);
            arma::vec lower_vals = z_col.elem(lower_idx);
            lower_candidate = std::max(lower_candidate, lower_vals.max());
          }
          arma::uvec upper_idx = arma::find(Y.col(i) == (level + 1));
          if (!upper_idx.is_empty()) {
            arma::vec z_col = z0.slice(0).col(i);
            arma::vec upper_vals = z_col.elem(upper_idx);
            upper_candidate = std::min(upper_candidate, upper_vals.min());
          }
          double finite_lower = clamp_probit_bound(lower_candidate);
          double finite_upper = clamp_probit_bound(upper_candidate);
          if (finite_lower >= finite_upper) {
            current_thresh(i, level) = 0.5 * (finite_lower + finite_upper);
          } else {
            current_thresh(i, level) = R::runif(finite_lower, finite_upper);
          }
        }
      }

      for(int j = 0; j < n; ++j){
        int category = static_cast<int>(Y.col(i)[j]);
        category = std::max(category, 1);
        category = std::min(category, K);
        double lower = current_thresh(i , category - 1);
        double upper = current_thresh(i , category);
        z0.slice(0).row(j).col(i) = rtnorm(mm(j), sd, lower, upper);
      }

      recenter_column_and_thresholds(z0, i, current_thresh);
      current_thresh(i, 0) = -arma::datum::inf;
      current_thresh(i, K) = arma::datum::inf;
    }

    for (int j = 0; j < k; ++j) {
      thresh.slice(j).row(s) = current_thresh.row(j);
    }

    for(int i = 0; i < k; ++i){
      D.row(i).col(i) = sqrt(1 / R::rgamma((delta + k - 1) / 2,
                   2 / arma::as_scalar(Rinv.slice(0).row(i).col(i))));
    }

    // expand latent data
    w = z0.slice(0) * D;

    // coefficients
    M = Sinv_X * X.t() * w;

    // gamma
    gamma = reshape(mvnrnd(reshape(M, k * p , 1),
                           kron(Sigma.slice(0), Sinv_X)),
                           p, k);

    // expanded coefficients
    beta  = gamma * Dinv.slice(0);

    // update Yhat
    Xbhat.slice(0) = X * beta;

    // error scatter matrix
    S_Y =   w.t() * w + I_k - M.t() * S_X * M;

    // } while (det(S_Y) < 0);

    // sample Psi
    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

    // sample Theta
    Theta.slice(0) =   wishrnd(inv(S_Y + Psi.slice(0)),  (deltaMP + k - 1) + (n - 1));

    // sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    // update Dinv
    Dinv.slice(0)  = inv(diagmat(sqrt(Sigma.slice(0).diag())));

    // inverse correlation matrix
    Rinv.slice(0)   = inv(cors);

    // update correlation matrix
    R.slice(0) = cors;

    beta_mcmc.slice(s) =reshape(beta, p,k);
    pcors_mcmc.slice(s) =  -(pcors - I_k);
    // cors_mcmc.slice(s) =  cors;
    // Sigma_mcmc.slice(s) = Sigma.slice(0);
    // Theta_mcmc.slice(s) = Theta.slice(0);
    // thresh.row(s) = thresh.slice(0).row(s);

  }

  arma::cube fisher_z = atanh(pcors_mcmc);

  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter - 50), 2);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] = pcor_mat;
  ret["beta"] = beta_mcmc;
  ret["thresh"]  = thresh;
  ret["fisher_z"] = fisher_z;
  return  ret;


}


// mixed data sampler
// [[Rcpp::export]]
Rcpp::List  copula(arma::mat z0_start,
                   arma::mat levels,
                   arma::vec K,
                   arma::mat Sigma_start,
                   int iter,
                   float delta,
                   float epsilon,
                   arma::vec idx,
                   bool progress
                   ) {

  // adapted from hoff 2008 for Bayesian hypothesis testing
  // with the matrix-F prior distribution for Theta

  // z0: latent data
  // levels: data matrix as sorted levels
  // K: levels in each columns

  Progress  pr(iter, progress);

  // number of rows
  float n = z0_start.n_rows;

  // number of columns
  int k = z0_start.n_cols;

  // k by k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  int nu = 1/ epsilon;
  // // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1;
  //
  // // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1;

  arma::uvec where ;

  // latent update
  arma::cube z0(n, k, 1,  arma::fill::zeros);

  z0.slice(0) = z0_start;
  //
  arma::cube zmcmc(n, k, iter, arma::fill::zeros);
  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  // scatter matrix dependent variables
  arma::mat S_Y(k, k, arma::fill::zeros);
  arma::mat B(epsilon * I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  Sigma.slice(0) = Sigma_start;

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  Theta.slice(0) = inv(Sigma_start);
  // arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);
//
  // covariance matrix
  // Sigma.slice(0) = Sigma_start;
  // arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  arma::vec lb(1);
  arma::vec ub(1);

  arma::mat mm(n,1);
  arma::mat ss(1,1);

  for(int  s = 1; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    for(int i = 0; i < k; ++i){

      mm = Sigma_i_not_i(Sigma.slice(0), i) *
        inv(remove_row(remove_col(Sigma.slice(0), i), i)) *
        remove_col(z0.slice(0), i).t();

      ss = select_row(Sigma.slice(0), i).col(i) -
        Sigma_i_not_i(Sigma.slice(0), i) *
        inv(remove_row(remove_col(Sigma.slice(0), i), i)) *
        Sigma_i_not_i(Sigma.slice(0), i).t();

      // sample latent data (0  = assumed continuous)
      if(idx(i) == 1){

        for(int r = 1; r  < K[i]+1; ++r){

          arma::uvec where = find(levels.col(i) == r );

          arma::mat temp1 = arma::conv_to<arma::mat>::from(where);

          int r_levels = temp1.n_elem;

          arma::vec lb_check =  {select_col(z0.slice(0), i).elem(find(levels.col(i) == (r - 1)))};
          arma::vec ub_check =  {select_col(z0.slice(0), i).elem(find(levels.col(i) == (r + 1)))};

          // Equation X in X
          if(lb_check.n_elem == 0){
            lb.fill(-arma::datum::inf);
          } else {
            lb.fill(lb_check.max());
          }

          if(ub_check.n_elem == 0){
            ub.fill(arma::datum::inf);
          } else {
            ub.fill(ub_check.min());
          }

          for(int l = 0; l < r_levels; ++l){

            z0.slice(0).col(i).row(temp1(l)) = R::qnorm(R::runif(
              R::pnorm(arma::as_scalar(lb),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE),
              R::pnorm(arma::as_scalar(ub),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE)),
              mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE);
          }
        }
      }
    }

    // novel matrix-F prior distribution
    // scatter matrix
    S_Y = z0.slice(0).t() * z0.slice(0);

    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

    // sample Theta
    Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

    // sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    pcors_mcmc.slice(s) =  -(pcors - I_k);
  }

  arma::cube fisher_z = atanh(pcors_mcmc);
  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter - 50), 2);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] = pcor_mat;
  ret["fisher_z"] = fisher_z;
  return ret;
}


// partials to correlations
// [[Rcpp::export]]
Rcpp::List pcor_to_cor_internal(arma::cube x, int p) {

  // slices
  int iter = x.n_slices;

  // correlation matrix (R)
  arma::cube R(p, p, iter);

  // precision matrix
  arma::mat Theta_s(p, p);

  // sigma matrix
  arma::mat Sigma_s(p, p);

  for(int s = 0; s < iter; ++s){

    Theta_s =   x.slice(s);

    for(int j = 0; j < p; ++j) {

      Theta_s.col(j).row(j) = 1;

    }

    arma::mat Sigma_s = inv(Theta_s);

    R.slice(s) =  diagmat(1 / sqrt(Sigma_s.diag())) * Sigma_s * diagmat(1 / sqrt(Sigma_s.diag()));

  }

  arma::mat R_mean = mean(R, 2);
  Rcpp::List ret;
  ret["R"] = R;
  ret["R_mean"] = R_mean;
  return ret;
}


// partials to correlations
// [[Rcpp::export]]
Rcpp::List predictability_helper(arma::mat Y,
                                 arma::colvec y,
                                 arma::cube XX,
                                 arma::mat Xy,
                                 int n,
                                 int iter) {


  arma::mat r2(iter, 1, arma::fill::zeros);
  arma::mat ppc(n,1);

  for(int s = 0; s < iter; ++s){

    // coefs
    arma::mat coefs = Xy.col(s).t() * inv(XX.slice(s));

    // predict mean
    arma::mat mu = Y * coefs.t();

    for(int j = 0; j < n; ++j){

      // ppc person j
      arma::vec ppc_j = Rcpp::rnorm(1, mu[j], stddev(mu - y));

      // y rep
      ppc.row(j) = arma::as_scalar(ppc_j);

    }

    // bayes R2
    r2.row(s) = var(mu) / var(ppc);

  }

  Rcpp::List ret;
  ret["r2"] = r2;
  return ret;
}

// regression coefficients
// [[Rcpp::export]]
Rcpp::List beta_helper_fast(arma::cube XX,
                       arma::mat Xy,
                       int p,
                       int iter) {


  arma::mat coefs(iter, p, arma::fill::zeros);

  for(int s = 0; s < iter; ++s){

    arma::vec coef_s = inv(XX.slice(s)).t() * Xy.col(s);
    coefs.row(s) =  coef_s.t();

  }

  Rcpp::List ret;
  ret["coefs"] = coefs;
  return ret;
}


// [[Rcpp::export]]
Rcpp::List pred_helper_latent(arma::mat Y,
                       arma::cube XX,
                       arma::mat Xy,
                       arma::vec quantiles,
                       int n,
                       int iter) {


  arma::mat yhat(iter, n, arma::fill::zeros);

  for(int s = 0; s < iter; ++s){

    arma::vec yhat_s = Y * inv(XX.slice(s)).t() * Xy.col(s) ;

    yhat.row(s) =  yhat_s.t();

  }

  // yhat
  arma::mat yhat_mean = mean(yhat);

  // quantiles
  arma::mat yhat_quantiles = quantile(yhat, quantiles);

  arma::mat yhat_sd = stddev(yhat);

  // returned
  Rcpp::List ret;
  ret["yhat"] = yhat;
  ret["yhat_mean"] = yhat_mean;
  ret["yhat_sd"] = yhat_sd;
  ret["yhat_quantiles"] = yhat_quantiles;
  return ret;
}


// [[Rcpp::export]]
float KL_univariate(float var_1, float var_2){

  float kl = log(sqrt(var_2)/sqrt(var_1)) + (var_1/(2 * var_2)) -  0.5;
  return kl;

}

// [[Rcpp::export]]
Rcpp::List ppc_helper_nodewise_fast(arma::cube Theta,
                                    int n1,
                                    int n2,
                                    int p){

  int iter = Theta.n_slices;

  Progress  pr(iter, TRUE);

  arma::vec mu(p, arma::fill::zeros);

  arma::mat kl(iter, p,  arma::fill::zeros);

  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    arma::mat Sigma = inv(Theta.slice(s));

    arma::mat R =  diagmat(1 / sqrt(Sigma.diag())) * Sigma * diagmat(1 / sqrt(Sigma.diag()));

    arma::mat Yrep_1  = mvnrnd(mu, R, n1).t();
    arma::mat Yrep_2  = mvnrnd(mu, R, n2).t();

    arma::mat Yrep_cor_1 = cor(Yrep_1);
    arma::mat Yrep_cor_2 = cor(Yrep_2);

    for(int j = 0; j < p; ++j){

      arma::mat pred_1 =  remove_col(Yrep_1, j) * trans(Sigma_i_not_i(Yrep_cor_1, j) *
        inv(remove_row(remove_col(Yrep_cor_1, j), j)));

      arma::mat pred_2 =  remove_col(Yrep_2, j) * trans(Sigma_i_not_i(Yrep_cor_2, j) *
        inv(remove_row(remove_col(Yrep_cor_2, j), j)));

      arma::mat var_1 = var(pred_1);
      arma::mat var_2 = var(pred_2);

      kl(s, j) = (KL_univariate(arma::as_scalar(var_1), arma::as_scalar(var_2)) +
        KL_univariate(arma::as_scalar(var_2), arma::as_scalar(var_1)))  * 0.5;

    }
  }

  Rcpp::List ret;
  ret["kl"] = kl;
  return ret;

}


// [[Rcpp::export]]
double KL_divergnece_mvn(arma::mat Theta_1, arma::mat Theta_2) {

  // number of variables
  int p = Theta_1.n_cols;

  arma::mat Sigma_1 = inv(Theta_1);

  // identity matrix
  arma::mat  I_p(p, p, arma::fill::eye);

  double kl = 0.50 * (trace(Sigma_1 * Theta_2) -  log(det(Sigma_1 * Theta_2)) - p);

  return kl;


}

// [[Rcpp::export]]
float sum_squares(arma::mat Rinv_1, arma::mat Rinv_2){

  arma::rowvec ss = sum(square(Rinv_1 - Rinv_2), 0);

  return sum(ss);

}

// [[Rcpp::export]]
arma::vec my_dnorm( arma::vec x, arma::vec means, arma::vec sds){

  int n = x.size();

  arma::vec res(n);

  for( int i=0; i < n; i++) res[i] = R::dnorm(x[i], means[i], sds[i],  FALSE) ;

  return res;
}

// [[Rcpp::export]]
float hamming_distance(arma::mat Rinv_1,
                       arma::mat Rinv_2,
                       float df1,
                       float df2,
                       float dens,
                       int pcors,
                       float BF_cut){

  // approximate post sd
  arma::mat se_1  = sqrt((1 - square(Rinv_1)) / (df1));
  arma::mat se_2  = sqrt((1 - square(Rinv_2)) / (df2));

  // upper-triangular
  arma::uvec ids = arma::trimatu_ind(size(se_1), 1);

  // partial correlations
  arma::vec r_1 = Rinv_1(ids);
  arma::vec r_2 = Rinv_2(ids);

  // sds
  arma::vec se_up_1 = se_1(ids);
  arma::vec se_up_2 = se_2(ids);

  // matrix of zeros
  arma::vec zerovec(pcors, arma::fill::zeros);

  // mat for 0's and 1's
  arma::vec sig_1(pcors, arma::fill::zeros);
  arma::vec sig_2(pcors, arma::fill::zeros);

  // density at zero
  arma::vec dens_1 = my_dnorm(zerovec, r_1, se_up_1);
  arma::vec dens_2 = my_dnorm(zerovec, r_2, se_up_2);

  //
  for(int i = 0; i < pcors; ++i){

    // check BF_cut (group 1)
    if((1 / (dens_1(i) / dens))  > BF_cut){

      sig_1(i) = 1;

    } else {

      sig_1(i) = 0;

    }

    // check BF_cut (group 2)
    if((1 / (dens_2(i) / dens))  > BF_cut){

      sig_2(i) = 1;

    } else {

      sig_2(i) = 0;

    }
  }

  return sum(square(sig_1 - sig_2));


}

// [[Rcpp::export]]
float correlation(arma::mat Rinv_1,
                  arma::mat Rinv_2){

  arma::uvec ids = arma::trimatu_ind(size(Rinv_1), 1);

  arma::vec r_1 = Rinv_1(ids);

  arma::vec r_2 = Rinv_2(ids);

  arma::mat cor_nets = cor(r_1, r_2);

  float cors = arma::as_scalar(cor_nets);

  return cors;

}


// [[Rcpp::export]]
Rcpp::List ppc_helper_fast(arma::cube Theta,
                           int n1,
                           int n2,
                           int p,
                           float BF_cut,
                           float dens,
                           bool ppc_ss,
                           bool ppc_cors,
                           bool ppc_hd){

  int iter = Theta.n_slices;

  // progress bar
  Progress  pr(iter, TRUE);

  // mean vector
  arma::vec mu(p, arma::fill::zeros);

  // KL storage
  arma::vec kl(iter, arma::fill::zeros);

  // sum of squares storage
  arma::vec ss(iter, arma::fill::zeros);

  // hamming distance storage
  arma::vec hd(iter, arma::fill::zeros);

  // correlation storage
  arma::vec cors(iter, arma::fill::zeros);

  arma::mat Yrep_Rinv_1(p, p);

  int df1 = n1 - (p-2) - 2;

  int df2 = n2 - (p-2) - 2;

  int pcors = (p * (p - 1)) * 0.5;

  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    arma::mat Sigma = inv(Theta.slice(s));

    arma::mat R =  diagmat(1 / sqrt(Sigma.diag())) * Sigma * diagmat(1 / sqrt(Sigma.diag()));

    arma::mat Yrep_1  = mvnrnd(mu, R, n1).t();
    arma::mat Yrep_2  = mvnrnd(mu, R, n2).t();

    arma::mat Yrep_Theta_1 = inv(cov(Yrep_1));
    arma::mat Yrep_Theta_2 = inv(cov(Yrep_2));

    arma::mat Yrep_Rinv_1 = diagmat(1 / sqrt(Yrep_Theta_1.diag())) *
      Yrep_Theta_1 *
      diagmat(1 / sqrt(Yrep_Theta_1.diag()));

    arma::mat Yrep_Rinv_2 = diagmat(1 / sqrt(Yrep_Theta_2.diag())) *
      Yrep_Theta_2 *
      diagmat(1 / sqrt(Yrep_Theta_2.diag()));

    kl(s) = 0.5 * (KL_divergnece_mvn(Yrep_Rinv_1, Yrep_Rinv_2) +
      KL_divergnece_mvn(Yrep_Rinv_2, Yrep_Rinv_1));

    if(ppc_ss){

      ss(s) = sum_squares(Yrep_Rinv_1, Yrep_Rinv_2) * 0.50;

    }

    if(ppc_hd){

      hd(s) = hamming_distance(Yrep_Rinv_1,
         Yrep_Rinv_2,
         df1, df2, dens,
         pcors, BF_cut);
    }

    if(ppc_cors){

      cors(s) = correlation(Yrep_Rinv_1, Yrep_Rinv_2);

    }

  }

  Rcpp::List ret;
  ret["kl"] = kl;
  ret["ss"] = ss;
  ret["hd"] = hd;
  ret["cors"] = cors;
  return ret;

}


// [[Rcpp::export]]
arma::mat mvnrnd(int n, arma::vec mu, arma::mat Sigma){

  arma::mat Y = mvnrnd(mu, Sigma, n).t();

  return Y;

}


// [[Rcpp::export]]
Rcpp::List var(arma::mat Y,
                         arma::mat X,
                         float delta,
                         float epsilon,
                         arma::mat beta_prior,
                         int iter,
                         arma::mat start,
                         bool progress){


  // progress
  Progress  pr(iter, progress);

  // number of rows
  int n = Y.n_rows;

  // number of dependent variables
  int k = Y.n_cols;

  // number of predictors
  int p = X.n_cols;

  int nu = 1/ epsilon;

  // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1;

  // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

  // k * k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  // p * p identity mat
  arma::mat  I_p(p, p, arma::fill::eye);

  // scatter matrix X' * X
  arma::mat S_X(X.t() * X + beta_prior);

  // inv S_X
  arma::mat Sinv_X(inv(S_X));

  // Psi update
  arma::cube Psi(k, k, 1, arma::fill::zeros);

  // scatter matrix dependent variables
  arma::mat S_Y(k, k, arma::fill::zeros);
  arma::mat B(epsilon*I_k);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  // precison matrix
  arma::cube Theta(k, k, 1, arma::fill::zeros);
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);

  // covariance matrix
  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  // coefficients
  arma::mat beta(p, k, arma::fill::zeros);
  arma::cube beta_mcmc(p, k, iter,  arma::fill::zeros);

  // starting value
  Sigma.slice(0) = inv(start);
  Psi.slice(0).fill(arma::fill::eye);
  Theta.slice(0) = start;

  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    // draw coefficients
    beta = reshape(mvnrnd(reshape(Sinv_X * X.t() * Y, k*p , 1),
                          kron(Sigma.slice(0), Sinv_X)),
                          p, k);

    // scatter matrix
    S_Y = Y.t() * Y + I_k - beta.t() * S_X * beta;

    // sample Psi
    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + k - 1);

    // sample Theta
    Theta.slice(0) =   wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + k - 1) + (n - 1));

    // Sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));


    beta_mcmc.slice(s) = beta;
    pcors_mcmc.slice(s) =  -(pcors - I_k);
  }

  arma::cube fisher_z = atanh(pcors_mcmc);

  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter - 50), 2);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] =  pcor_mat;
  ret["beta"] = beta_mcmc;
  ret["fisher_z"] = fisher_z;
  return ret;
}


// [[Rcpp::export]]
Rcpp::List hft_algorithm(arma::mat Sigma, arma::mat adj, double tol, double max_iter) {

  arma::mat S = Sigma;
  arma::mat W = S;

  arma::uvec upper_indices = trimatu_ind( size(S) );
  arma::mat W_previous = S;
  double p = S.n_cols;
  arma::mat iter(1,1, arma::fill::zeros);
  double max_diff = 100;
  arma::mat w_12(1, p-1);

  while(max_diff > tol){

    for(int i = 0; i < p; ++i){

      arma::mat beta(1,p-1, arma::fill::zeros);
      arma::uvec pad_index =  find(Sigma_i_not_i(adj, i) == 1);

      if(pad_index.n_elem == 0 ){
        w_12 = beta;
      } else {

        arma::mat W_11 = remove_row(remove_col(W , i), i);
        arma::mat s_12 = Sigma_i_not_i(S,i);

        arma::mat W_11_star = W_11.submat(pad_index, pad_index);
        arma::mat s_12_star = s_12(pad_index);

        beta(pad_index) = inv(W_11_star) * s_12_star;
        arma::mat w_12 = W_11 * beta.t();
        arma::mat temp = W.col(i).row(i);
        w_12.insert_rows(i, temp);

        for(int k = 0; k < p; ++k){
          W.row(i).col(k) = w_12(k);
          W.row(k).col(i) = w_12(k);
        }

        max_diff = max(W.elem(upper_indices) -  W_previous.elem(upper_indices));
        W_previous = W;
      }
    }

    iter(0,0) = iter(0,0) + 1;

    if(iter(0,0) == max_iter){
      break;
    }

  }

  arma::mat Theta = inv(W) % adj;
  Rcpp::List ret;
  ret["Theta"] = Theta;
  return ret;
}

// [[Rcpp::export]]
double bic_fast(arma::mat Theta, arma::mat S, double n, float prior_prob){

  arma::mat UU = trimatu(Theta,  1);

  arma::vec nonzero = nonzeros(UU);
  double neg_ll =  -2 * ((n*0.5) * (log(det(Theta)) - trace(S * Theta)));

  double bic = neg_ll + (nonzero.n_elem * log(n) - (nonzero.n_elem * log(prior_prob / (1 - prior_prob))));

  return bic;
}

// [[Rcpp::export]]
Rcpp::List find_ids(arma::mat x){
  arma::mat UU = trimatu(x,  1);
  arma::uvec alt_lower_indices = trimatl_ind( size(x),  -1);
  UU.elem(alt_lower_indices).fill(10);
  UU.diag().fill(10);
  arma::uvec zero = find(UU == 0);
  arma::uvec nonzero = find(UU == 1);

  Rcpp::List ret;
  ret["nonzero"] = nonzero;
  ret["zero"] = zero;
  return ret;
}

// Search for possible graphs, accepting if BIC_new is better than BIC_old
// MH algo does not work probabilistically - accept only better, never worse.
// Revisit this function 
// [[Rcpp::export]]
Rcpp::List search(arma::mat S,
                  float iter,
                  double old_bic,
                  arma::mat start_adj,
                  float n,
                  float gamma,
                  int stop_early,
                  bool progress){

  Progress  pr(iter, progress);

  int p = S.n_cols;

  arma::cube adj(p, p, iter);
  
  // Copy start_adj to adj_s
  arma::mat adj_s = start_adj; 

  // Create start object containing position of zero and nonzeros
  Rcpp::List start = find_ids(start_adj);
  
  arma::uvec zeros = start["zero"];
  arma::uvec nonzeros = start["nonzero"];

  // Use adj_start from R
  arma::mat mat_old = start_adj;
  // Initialize adj_mat to match mat_old
  arma::mat adj_mat = mat_old;    
  
  // initialize vectors
  arma::vec bics(iter, arma::fill::zeros);
  arma::vec acc(1, arma::fill::zeros);
  arma::vec repeats(1, arma::fill::zeros);
 
  // Loop through iterations
  for(int s = 0; s < iter; ++s){

    // Incrementing progress bar
    pr.increment();

    // Catch user abort key
    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    adj_s = mat_old;

    // Flip one edge at a time
    if (s % 2 == 0){
      arma::vec id_add = Rcpp::RcppArmadillo::sample(arma::conv_to<arma::vec>::from(zeros), 1, false);
      adj_s.elem(arma::conv_to<arma::uvec>::from(id_add)).fill(1);
    } else {
      arma::vec id_add = Rcpp::RcppArmadillo::sample(arma::conv_to<arma::vec>::from(nonzeros), 1, false);
      adj_s.elem(arma::conv_to<arma::uvec>::from(id_add)).fill(0);
    }

    // Ensure that adj_mat is symmetric
    adj_mat = symmatu(adj_s);
    adj_mat.diag().fill(1);
    
    // Run the hft_algorithm and compute the BIC
    Rcpp::List fit1 = hft_algorithm(S, adj_mat, 0.00001, 10);
    arma::mat  Theta = fit1["Theta"];
    double new_bic = bic_fast(Theta, S, n, gamma);
    
    // Specifically compute delta to facilitate debugging
    double delta =  new_bic - old_bic;

    // Generate a random uniform number for probabilistic acceptance
    // double random_uniform = arma::randu();

    // Metropolis-Hastings acceptance criterion
    if(exp(-0.5 *  delta ) >= 1) { //random_uniform ){
      // go back to >=1, as accepting probabilistically has BIC creep up for unknown reason
      mat_old = adj_mat;
      adj.slice(s) = adj_mat;
      old_bic = new_bic;
      acc(0)++;
      repeats(0) = 0;
      Rcpp::List start =  find_ids(adj_mat);
      arma::uvec zeros = start["zero"];
      arma::uvec nonzeros = start["nonzero"];
    } else {
      adj.slice(s) = mat_old;
      repeats(0)++;
    }
    
    bics(s) = old_bic;

    if(repeats(0) > stop_early){
      Rcpp::Rcout << "Stopping early at iteration " << s << std::endl;
      break;
    } 
    
  }

  Rcpp::List ret;
  ret["p"] = p;
  ret["adj_mat"] = adj_mat;
  ret["bics"] = bics;
  ret["adj"]= adj;
  ret["acc"] = acc;

  return ret;
}


// random walk sampler for known graph using matrix-F

// [[Rcpp::export]]
Rcpp::List fast_g_matrix_F(arma::mat Y,
                           arma::mat adj,
                           arma::vec mu_samples,
                           arma::mat cov_samples,
                           int iter,
                           int p,
                           float N,
                           float prior_sd,
                           float kappa1,
                           bool progress){


  Progress  pr(iter, progress);

  arma::cube Theta_G(p, p, iter, arma::fill::zeros);

  arma::vec uniform(1, 1);

  arma::vec acc(1, arma::fill::zeros);

  arma::mat UU = trimatl(adj,  1);

  arma::uvec alt_lower_indices = trimatu_ind(size(adj),  1);

  UU.elem(alt_lower_indices).fill(10);

  UU.diag().fill(1);

  arma::uvec nonzero = find(UU == 1);

  arma::vec kappa_store(iter, arma::fill::zeros);

  arma::mat Theta_can1(p, p, arma::fill::zeros);

  arma::mat Theta_s1(p,p,arma::fill::zeros);

  Theta_s1.elem(nonzero) = mu_samples;

  arma::mat Theta_s = symmatl(Theta_s1);

  arma::mat SSY = Y.t() * Y;

  float log_det_s = log(det(Theta_s));

  arma::mat Binv(p, p, arma::fill::zeros);

  Binv.diag().fill(10000);

  arma::vec acc_prob(1, 1, arma::fill::zeros);

  arma::mat I_p(p, p, arma::fill::eye);

  float deltaF   = 1/(prior_sd * prior_sd) - 1;

  double logpriorF_s = (10000 - p - 1) / 2 *
                       log_det_s - (deltaF + 10000 + p - 1)/2 *
                       log(det(I_p + Theta_s * Binv));

  double log_lik_s = (N/2) * log_det_s - 0.5 * trace(SSY * Theta_s);

  double log_post_s = logpriorF_s + log_lik_s;

  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }


    arma::vec theta_can =  mvnrnd(mu_samples, kappa1 * cov_samples);

    Theta_can1.elem(nonzero) = theta_can.t();

    arma::mat Theta_can = symmatl(Theta_can1);

    if(Theta_can.is_sympd()) {

      float log_det_can =  log(det(Theta_can));

      double logpriorF_can = (10000 - p - 1) / 2 *
                             log_det_can - (deltaF + 10000 + p - 1) / 2 *
                             log(det(I_p + Theta_can * Binv));

      double log_lik_can = (N/2) * log_det_can - 0.5 * trace(SSY * Theta_can);

      double log_post_can = logpriorF_can + log_lik_can;

      arma::vec uniform(1, 1, arma::fill::randu);

      double test = exp(log_post_can -  log_post_s);

      if(test >  uniform(0) ){

        acc(0) = acc(0) + 1;

        Theta_s = Theta_can;

        mu_samples = theta_can;

        log_post_s = log_post_can;

      }

      acc_prob(0) = acc(0) / s;

      if(acc_prob(0) < 0.30){
        kappa1 = kappa1 * 0.9;
      }

      if(acc_prob(0) > 0.50){
        kappa1 = kappa1 * 1.1;
      }
    }

    kappa_store(s) = kappa1;

    Theta_G.slice(s) = Theta_s;

  }

  Rcpp::List ret;
  ret["acc"] = acc;
  ret["Theta_G"] = Theta_G;
  ret["acc_prob"] = acc_prob;
  ret["kappa"] = kappa_store;
  return ret;
}


// [[Rcpp::export]]
arma::cube contrained_helper(arma::cube cors,
                             arma::mat adj,
                             int iter,
                             bool progress){


  Progress  pr(iter, progress);

  int p = cors.slice(0).n_cols;

  arma::cube Theta(p, p, iter, arma::fill::zeros);

  for(int s = 0; s < iter; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    Rcpp::List fit1 = hft_algorithm(cors.slice(s), adj, 0.00001, 10);

    arma::mat  Theta_s = fit1["Theta"];

    Theta.slice(s) = Theta_s;

    }

  return Theta;
}

// [[Rcpp::export]]
Rcpp::List missing_copula(arma::mat Y,
                             arma::mat Y_missing,
                             arma::mat z0_start,
                             arma::mat Sigma_start,
                             arma::mat levels,
                             int iter_missing,
                             bool progress_impute,
                             arma::vec K,
                             arma::vec idx,
                             float epsilon,
                             float delta) {
  // progress
  Progress  pr(iter_missing, progress_impute);

  int p = Y.n_cols;
  int n = Y.n_rows;

  arma::mat Y_impute = Y;

  arma::cube Y_collect(n, p, iter_missing, arma::fill::zeros);

  // p by p identity mat
  arma::mat  I_p(p, p, arma::fill::eye);
  arma::cube Psi(p, p, 1, arma::fill::zeros);

  // latent update
  arma::cube  z0(n, p, 1,  arma::fill::zeros);

  z0.slice(0) = z0_start;

  arma::uvec index = find(Y_missing == 1);

  //No more use:  int n_na = index.n_elem;

  int nu = 1/ epsilon;

  // // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + p - 1;

  // // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - p + 1;

  arma::mat B(epsilon * I_p);
  arma::mat BMP(inv(B));
  arma::mat BMPinv(inv(BMP));

  arma::mat z(n,p);

  arma::vec lb(1);
  arma::vec ub(1);

  arma::cube Sigma(p, p, 1, arma::fill::zeros);
  arma::cube Theta(p, p, 1, arma::fill::zeros);
  Sigma.slice(0) = Sigma_start;

  arma::mat mm(n,1);
  arma::mat ss(1,1);

  // partial correlations
  arma::mat pcors(p,p);
  arma::cube pcors_mcmc(p, p, iter_missing, arma::fill::zeros);

  for(int  s = 0; s < iter_missing; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    for(int i = 0; i < p; ++i){

      mm = Sigma_i_not_i(Sigma.slice(0), i) *

        inv(remove_row(remove_col(Sigma.slice(0), i), i)) *

        remove_col(z0.slice(0), i).t();

      ss = select_row(Sigma.slice(0), i).col(i) -
        Sigma_i_not_i(Sigma.slice(0), i) *
        inv(remove_row(remove_col(Sigma.slice(0), i), i)) *
        Sigma_i_not_i(Sigma.slice(0), i).t();

      // sample latent data (0  = assumed continuous)
      if(idx(i) == 1){

        for(int r = 1; r  < K[i]+1; ++r){

          arma::uvec where = find(levels.col(i) == r);

          arma::mat temp1 = arma::conv_to<arma::mat>::from(where);

          int r_levels = temp1.n_elem;

          arma::vec lb_check =  {select_col(z0.slice(0), i).elem(find(levels.col(i) == (r - 1)))};
          arma::vec ub_check =  {select_col(z0.slice(0), i).elem(find(levels.col(i) == (r + 1)))};

          // sample accoring to ranks
          if(lb_check.n_elem == 0){
            lb.fill(-arma::datum::inf);
          } else {
            lb.fill(lb_check.max());
          }

          if(ub_check.n_elem == 0){
            ub.fill(arma::datum::inf);
          } else {
            ub.fill(ub_check.min());
          }

          for(int l = 0; l < r_levels; ++l){

            z0.slice(0).col(i).row(temp1(l)) = R::qnorm(R::runif(
              R::pnorm(arma::as_scalar(lb),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE),
              R::pnorm(arma::as_scalar(ub),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE)),
              mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE);
          }

        }

      }

      arma::vec Y_j = Y_missing.col(i);

      double check_na = sum(Y_j);

      if(check_na > 0){

        arma::uvec  index_j = find(Y_missing.col(i) == 1);

        int  n_missing = index_j.n_elem;

        for(int m = 0; m < n_missing;  ++m){

          arma::vec ppd_i = Rcpp::rnorm(1,  mm(index_j[m]), sqrt(ss(0)));

          z0.slice(0).col(i).row(index_j[m]) = ppd_i(0);

        }

      }

    }

    arma::mat S_Y = z0.slice(0).t() * z0.slice(0);

    Psi.slice(0) = wishrnd(inv(BMPinv + Theta.slice(0)), nuMP + deltaMP + p - 1);

    // sample Theta
    Theta.slice(0) = wishrnd(inv(Psi.slice(0) + S_Y),  (deltaMP + p - 1) + (n - 1));

    // Sigma
    Sigma.slice(0) = inv(Theta.slice(0));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    // store posterior samples
    pcors_mcmc.slice(s) =  -(pcors - I_p);

  }

  arma::cube fisher_z = atanh(pcors_mcmc);
  arma::mat  pcor_mat = mean(pcors_mcmc.tail_slices(iter_missing - 50), 2);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["pcor_mat"] = pcor_mat;
  ret["fisher_z"] = fisher_z;
  ret["Y_collect"] = Y_collect;
  return  ret;
}



// [[Rcpp::export]]
Rcpp::List missing_copula_data(arma::mat Y,
                   arma::mat Y_missing,
                   arma::mat z0_start,
                   arma::mat Sigma_start,
                   arma::mat levels,
                   int iter_missing,
                   bool progress_impute,
                   arma::vec K,
                   arma::vec idx,
                   float lambda) {
  // progress
  Progress  pr(iter_missing, progress_impute);

  int p = Y.n_cols;

  int n = Y.n_rows;

  arma::mat Y_impute = Y;

  arma::cube Y_collect(n, p,
                       iter_missing,
                       arma::fill::zeros);

  arma::mat I_p(p, p, arma::fill::eye);

  arma::cube z0(n, p, 1,  arma::fill::zeros);

  z0.slice(0) = z0_start;

  arma::vec lb(1);
  arma::vec ub(1);

  arma::cube Sigma(p, p, 1, arma::fill::zeros);
  Sigma.slice(0) = Sigma_start;

  arma::mat mm(n,1);
  arma::mat ss(1,1);

  for(int  s = 0; s < iter_missing; ++s){

    pr.increment();

    if (s % 250 == 0){
      Rcpp::checkUserInterrupt();
    }

    for(int i = 0; i < p; ++i){

      mm = Sigma_i_not_i(Sigma.slice(0), i) *
        inv(remove_row(remove_col(Sigma.slice(0), i), i)) *
        remove_col(z0.slice(0), i).t();

      ss = select_row(Sigma.slice(0), i).col(i) -
        Sigma_i_not_i(Sigma.slice(0), i) *
        inv(remove_row(remove_col(Sigma.slice(0), i), i)) *
        Sigma_i_not_i(Sigma.slice(0), i).t();

      if(idx(i) == 1){

        for(int r = 1; r  < K[i] + 1; ++r){

          arma::uvec where = find(levels.col(i) == r);

          arma::mat temp1 = arma::conv_to<arma::mat>::from(where);

          int r_levels = temp1.n_elem;

          arma::vec lb_check =  {select_col(z0.slice(0), i).elem(find(levels.col(i) == (r - 1)))};
          arma::vec ub_check =  {select_col(z0.slice(0), i).elem(find(levels.col(i) == (r + 1)))};

          if(lb_check.n_elem == 0){

            lb.fill(-arma::datum::inf);

          } else {

            lb.fill(lb_check.max());

          }

          if(ub_check.n_elem == 0){

            ub.fill(arma::datum::inf);

          } else {

            ub.fill(ub_check.min());

          }

          for(int l = 0; l < r_levels; ++l){

            z0.slice(0).col(i).row(temp1(l)) = R::qnorm(R::runif(
              R::pnorm(arma::as_scalar(lb),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE),
              R::pnorm(arma::as_scalar(ub),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE)),
              mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE);
          }
        }
      }

      arma::vec Y_j = Y_missing.col(i);

      double check_na = sum(Y_j);


      if(check_na > 0){

        arma::uvec  index_j = find(Y_missing.col(i) == 1);

        int  n_missing = index_j.n_elem;

        for(int m = 0; m < n_missing; ++m){

          arma::vec ppd_i = Rcpp::rnorm(1,  mm(index_j[m]), sqrt(ss(0)));

          z0.slice(0).col(i).row(index_j[m]) = ppd_i(0);

        }
      }
    }

    arma::mat S_Y = z0.slice(0).t() * z0.slice(0);

    Sigma.slice(0) =  inv(wishrnd(inv(S_Y + I_p * lambda), n + lambda));

    for(int i = 0; i < p; ++i){

      arma::mat  zz = z0.slice(0).col(i);

      ss = Sigma.slice(0).col(i).row(i);

      arma::vec Y_jnew = Y_missing.col(i);

      double check_nanew = sum(Y_jnew);

      if(check_nanew > 0){

        arma::uvec  index_j = find(Y_missing.col(i) == 1);

        int  n_missing = index_j.n_elem;

        for(int m = 0; m < n_missing; ++m){

          Y_impute.col(i).row(index_j[m]) =  quantile_type_1(Y.col(i),
                       R::pnorm(zz(index_j[m]), 0, sqrt(ss(0)), TRUE, FALSE));
          }
      }
    }

    Y_collect.slice(s) = Y_impute;

  }

  Rcpp::List ret;
  ret["Y_collect"] = Y_collect;
  return  ret;
}
