// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <progress.hpp>
#include <progress_bar.hpp>
#include <truncnorm.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppProgress)]]

// helpers are first (avoids separate files)


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
        arma::vec ppc_i = Rcpp::rnorm(1,  pred(index_j[i]), arma::conv_to<double>::from(sd_j));
        Y.col(j).row(index_j[i]) = arma::conv_to<double>::from(ppc_i);
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
                            bool store_all) {


  // progress
  Progress  pr(iter_missing, progress_impute);

  int p = Y.n_cols;
  int n = Y.n_rows;

  arma::uvec index = find(Y_missing == 1);

  int n_na = index.n_elem;

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
                                      arma::conv_to<double>::from(sd_j));

        Y.col(j).row(index_j[i]) = arma::conv_to<double>::from(ppd_i);

      }

    }

    arma::mat S_Y = Y.t() * Y;
    arma::mat Theta = wishrnd(inv(S_Y), (n - 1));
    Sigma = inv(Theta);
    ppd_missing.row(s) = Y.elem(index).t();

    if(store_all){
      Y_all.slice(s) = Y;
    }
  }

  arma::vec lb = {0.025};
  arma::vec ub = {0.975};

  arma::mat  ppd_mean = mean(ppd_missing, 0).t();
  arma::mat  ppd_sd = stddev(ppd_missing, 0, 0).t();
  arma::mat  ppd_lb = quantile(ppd_missing, lb).t();
  arma::mat  ppd_ub = quantile(ppd_missing, ub).t();

  arma::mat  ppd_summary = join_rows(ppd_mean, ppd_sd, ppd_lb, ppd_ub);

  Rcpp::List ret;
  ret["Y"] = Y;
  ret["Y_all"] = Y_all;
  ret["ppd_missing"] = ppd_missing;
  ret["ppd_mean"] = ppd_mean;
  ret["ppd_summary"] = ppd_summary;
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
                                               iter_missing = iter_missing);
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

    a2.row(i) = a.col(arma::conv_to<int>::from(y.row(i)));

    b2.row(i) = b.col(arma::conv_to<int>::from(y.row(i)));

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
                   2 / arma::conv_to<double>::from(Rinv.slice(0).row(i).col(i))));
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


    Dinv.slice(0)  = inv(diagmat(sqrt(Sigma.slice(0).diag())));

    Rinv.slice(0)   = inv(cors);

    R.slice(0) = cors;


    // // latent data
    // for(int i = 0; i < n; ++i){
    //
    //   Rcpp::List z_samples = trunc_mvn(Xbhat.slice(0).row(i).t(),
    //                                    Rinv.slice(0),
    //                                    z0.slice(0).row(i).t(),
    //                                    Y.row(i).t(),
    //                                    cutpoints);
    //
    //   arma::mat z_i = z_samples[0];
    //
    //   z0.slice(0).row(i) = z_i.t();
    //
    // }

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

  // story thresholds
  arma::cube thresh_mcmc(k, n_levels.size() + 1, iter);

  // MH R ratio
  arma::mat R_numer(n, k);
  arma::mat R_denom(n, k);

  arma::mat mm(n,1);
  arma::mat ss(1,1);


  // draw coefs conditional on w
  arma::mat gamma(p, k, arma::fill::zeros);

  // starting thresholds
  for(int i = 0; i < max(Y.col(0)) - 1; ++i){
    for(int j = 0; j < k; ++j){
      thresh_mat(j, i + 2, 0) = R::qnorm5(sum(Y.col(j) <= thresh_sampled[i]) / n,
                 -R::qnorm5(sum(Y.col(j) == 1) / n, 0, 1, TRUE, FALSE), 1, TRUE, FALSE);
    }
  }

  // start sampling
  for(int s = 0; s < iter; ++s){

    Rcpp::checkUserInterrupt();

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
      for(int j = 0; j < n; ++j){

        z0.slice(0).row(j).col(i) = R::qnorm(R::runif(
          // minimum
          R::pnorm(thresh_mat.slice(0)(i , (Y.col(i)[j] - 1)), mm(j), sqrt(ss(0)), TRUE, FALSE),

          // maximum
          R::pnorm(thresh_mat.slice(0)(i , Y.col(i)[j]), mm(j), sqrt(ss(0)), TRUE, FALSE)),

          // location and scale
          mm(j), sqrt(ss(0)), TRUE, FALSE);
      }
    }

    for(int i = 0; i < max(Y.col(0)) - 1; ++i){
      for(int j = 0; j < k ; ++j){
        arma::mat slice_c_thresh = thresh_mat.slice(0);
        c_thresh_mat(j, i+2, 0) = R::qnorm5(R::runif(R::pnorm(0,  slice_c_thresh(j , thresh_sampled[i]), MH, TRUE, FALSE), 1),
                     slice_c_thresh(j , thresh_sampled[i]), MH, TRUE, FALSE);
      }
    }

    for(int j = 0; j < n; ++j){
      for(int i = 0; i < k; ++i){

        R_numer(j, i) =  (R::pnorm5(c_thresh_mat(i , Y.col(i)[j], 0) - mm(j), 0, 1, TRUE, FALSE) -
          R::pnorm5(c_thresh_mat(i , Y.col(i)[j]-1 , 0) - mm(j), 0, 1, TRUE, FALSE));

        R_denom(j, i) = (R::pnorm5(thresh_mat(i , Y.col(i)[j], 0) - mm(j), 0, 1, TRUE, FALSE) -
          R::pnorm5(thresh_mat(i , Y.col(i)[j] -1, 0) - mm(j), 0, 1, TRUE, FALSE));
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
        thresh_mat.slice(0).row(j) = c_thresh_mat.slice(0).row(j);
      }
    }

    for(int i = 0; i < k; ++i){
      D.row(i).col(i) = sqrt(1 / R::rgamma((deltaMP + k - 1) / 2,
                   2 / arma::conv_to<double>::from(Rinv.slice(0).row(i).col(i))));
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
    thresh_mcmc.slice(s) = thresh_mat.slice(0);

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

// ordinal sampler
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

  for(int i = 0; i < k; ++i){
    thresh.slice(i).col(0).fill(-arma::datum::inf);
    thresh.slice(i).col(K).fill(arma::datum::inf);
  }

  for(int i = 0; i < (K-2); ++i){
    for(int j = 0; j < k; ++j){
      thresh.slice(j).col(i+2).fill(i+1);
    }
  }

  // store thresholds
  arma::mat thresh_mcmc(K+1,iter, arma::fill::zeros);

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

      if(s == 1){
        // generate latent data
        for(int j = 0; j < n; ++j){

          z0.slice(0).row(j).col(i) = R::qnorm(R::runif(
            // minimum
            R::pnorm(thresh.slice(i)(0 , (Y.col(i)[j] - 1)), mm(j), sqrt(ss(0)), TRUE, FALSE),

            // maximum
            R::pnorm(thresh.slice(i)(0 , Y.col(i)[j]), mm(j), sqrt(ss(0)), TRUE, FALSE)),

            // location and scale
            mm(j), sqrt(ss(0)), TRUE, FALSE);
        }


      } else{


        for(int i = 0; i < k; ++i){

          for(int j = 2; j < (K); ++j){
            arma::vec lb = {select_col(z0.slice(0), i).elem(find(Y.col(i) == j)).max(), thresh.slice(i)(s-1, j-1) };
            arma::vec ub = {select_col(z0.slice(0), i).elem(find(Y.col(i) == j+1)).min(), thresh.slice(i)(s-1, j+1) };
            arma::vec v = Rcpp::runif(1,  lb.max(), ub.min());
            thresh.slice(i).row(s).col(j) =  arma::conv_to<double>::from(v);

          }
        }

        // generate latent data
        for(int j = 0; j < n; ++j){

          z0.slice(0).row(j).col(i) = R::qnorm(R::runif(
            // minimum
            R::pnorm(thresh.slice(i)(s , (Y.col(i)[j] - 1)), mm(j), sqrt(ss(0)), TRUE, FALSE),

            // maximum
            R::pnorm(thresh.slice(i)(s  , Y.col(i)[j]), mm(j), sqrt(ss(0)), TRUE, FALSE)),

            // location and scale
            mm(j), sqrt(ss(0)), TRUE, FALSE);
        }
        }
      }

    for(int i = 0; i < k; ++i){
      D.row(i).col(i) = sqrt(1 / R::rgamma((delta + k - 1) / 2,
                   2 / arma::conv_to<double>::from(Rinv.slice(0).row(i).col(i))));
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
              R::pnorm(arma::conv_to<double>::from(lb),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE),
              R::pnorm(arma::conv_to<double>::from(ub),  mm(temp1(l)), sqrt(ss(0)), TRUE, FALSE)),
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

    // // correlation
    // cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
    //   Sigma.slice(0) *
    //   diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    pcors_mcmc.slice(s) =  -(pcors - I_k);
    // cors_mcmc.slice(s) =  cors;
    // Sigma_mcmc.slice(s) = Sigma.slice(0);
    // Theta_mcmc.slice(s) = Theta.slice(0);
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
      ppc.row(j) = arma::conv_to<double>::from(ppc_j);

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

      kl(s, j) = (KL_univariate(arma::conv_to<float>::from(var_1), arma::conv_to<float>::from(var_2)) +
        KL_univariate(arma::conv_to<float>::from(var_2), arma::conv_to<float>::from(var_1)))  * 0.5;

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

  float cors = arma::conv_to<float>::from(cor_nets);

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

