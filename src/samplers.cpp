// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <truncnorm.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


// matrix F continous sampler
// Williams, D. R., & Mulder, J. (2019). Bayesian hypothesis testing for Gaussian
// graphical models:  Conditional independence and order constraints.

// matrix F for the precision matrix
//' @title Testing mat F
//' @name Theta_continuous

// [[Rcpp::export]]
Rcpp::List Theta_continuous(arma::mat Y,
                            int iter,
                            float delta,
                            float epsilon,
                            int prior_only, int explore) {


  // note p changed to k to be consistent
  //with the multivariate regression samplers

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

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    // store posterior samples
    pcors_mcmc.slice(s) =  -(pcors - I_k);
    cors_mcmc.slice(s) =  cors;
    Sigma_mcmc.slice(s) = Sigma.slice(0);
    Theta_mcmc.slice(s) = Theta.slice(0);



  }

  arma::cube fisher_z = atanh(pcors_mcmc);

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["cors"] =  cors_mcmc;
  ret["Theta"] = Theta_mcmc;
  ret["Sigma"] = Sigma_mcmc;
  ret["fisher_z"] = fisher_z;
  return ret;
}





//' @title multivariate regression
//' @name mv_continuous
//' @export

// [[Rcpp::export]]
Rcpp::List mv_continuous(arma::mat Y,
                          arma::mat X,
                          float delta,
                          float epsilon,
                          int iter){

  // number of rows
  int n = Y.n_rows;

  // number of dependent variables
  int k = Y.n_cols;

  // number of predictors
  int p = X.n_cols;

  int nu = 1/ epsilon;
  // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;

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
  Sigma.slice(0).fill(arma::fill::eye);
  Psi.slice(0).fill(arma::fill::eye);
  Theta.slice(0).fill(arma::fill::eye);

  for(int s = 0; s < iter; ++s){

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

    // store posterior samples
    beta_mcmc.slice(s) = beta;
    pcors_mcmc.slice(s) =  -(pcors - I_k);
    cors_mcmc.slice(s) =  cors;
    Sigma_mcmc.slice(s) = Sigma.slice(0);
    Theta_mcmc.slice(s) = Theta.slice(0);
  }

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["cors"] =  cors_mcmc;
  ret["beta"] = beta_mcmc;
  ret["Theta"] = Theta_mcmc;
  ret["Sigma"] = Sigma_mcmc;
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


//' @title testing binary
//' @name mv_binary
//' @export

// binary sampler
// [[Rcpp::export]]
Rcpp::List mv_binary(arma::mat Y,
                      arma::mat X,
                      float delta,
                      float epsilon,
                      int iter,
                      float beta_prior,
                      arma::rowvec cutpoints){

  // Y: data matrix (n * k)
  // X: predictors (n * p) (blank for "network")
  // iter: number of samples
  // cutpoints: truncation points
  // delta: hyperparameter
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
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);

  // covariance matrix
  arma::cube Sigma(k, k, 1, arma::fill::zeros);
  arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  // latent data
  arma::cube z0(n, k, 1,  arma::fill::zeros);

  // expanded latent data
  arma::mat w(n, k, arma::fill::zeros);

  // conditonal data
  arma::cube Xbhat(n, k,  1, arma::fill::zeros);

  // Rinv update
  arma::cube Rinv(k, k, 1, arma::fill::zeros);

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
  Sigma.slice(0).fill(arma::fill::eye);
  Psi.slice(0).fill(arma::fill::eye);
  Theta.slice(0).fill(arma::fill::eye);
  Dinv.slice(0).fill(arma::fill::eye);
  Rinv.slice(0).fill(arma::fill::eye);


  // start sampling
  for(int s = 0; s < iter; ++s){

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

    // latent data
    for(int i = 0; i < n; ++i){

      Rcpp::List z_samples = trunc_mvn(Xbhat.slice(0).row(i).t(),
                                       Rinv.slice(0),
                                       z0.slice(0).row(i).t(),
                                       Y.row(i).t(),
                                       cutpoints);

      arma::mat z_i = z_samples[0];

      z0.slice(0).row(i) = z_i.t();

    }

    beta_mcmc.slice(s) =reshape(beta, p,k);
    pcors_mcmc.slice(s) =  -(pcors - I_k);
    cors_mcmc.slice(s) =  cors;
    Sigma_mcmc.slice(s) = Sigma.slice(0);
    Theta_mcmc.slice(s) = Theta.slice(0);
  }

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["cors"] =  cors_mcmc;
  ret["beta"] = beta_mcmc;
  ret["Theta"] = Theta_mcmc;
  ret["Sigma"] = Sigma_mcmc;
  return  ret;
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


//' @title testing ordinal cowles for thresholds
//' @name mv_ordinal_cowles
//' @export

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



//' @title testing ordinal albert and chibbs for thresholds
//' @name mv_ordinal_albert
//' @export

// ordinal sampler
// [[Rcpp::export]]
Rcpp::List mv_ordinal_albert(arma::mat Y,
                          arma::mat X,
                          int iter,
                          float delta,
                          float epsilon,
                          int K){


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

  Rinv.slice(0).fill(arma::fill::eye);
  R.slice(0).fill(arma::fill::eye);
  Dinv.slice(0).fill(arma::fill::eye);
  Psi.slice(0).fill(arma::fill::eye);
  Sigma.slice(0).fill(arma::fill::eye);

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

    Rcpp::checkUserInterrupt();

    for(int i = 0; i < k; ++i){

      mm = Xbhat.slice(0).col(i).t() +
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        (remove_col(z0.slice(0), i).t() - remove_col(Xbhat.slice(0), i).t());

      ss = select_row(R.slice(0), i).col(i) -
        Sigma_i_not_i(R.slice(0), i) *
        inv(remove_row(remove_col(R.slice(0), i), i)) *
        Sigma_i_not_i(R.slice(0), i).t();



      // Rcpp::checkUserInterrupt();

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
    Theta.slice(0) =   wishrnd(inv(  S_Y),  (deltaMP + k - 1) + (n - 1));

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
    // thresh.row(s) = thresh.slice(0).row(s);

  }

  ////////////////////////////////////////////

  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["cors"] =  cors_mcmc;
  ret["beta"] = beta_mcmc;
  ret["Theta"] = Theta_mcmc;
  ret["Sigma"] = Sigma_mcmc;
  ret["thresh"]  = thresh;
  return  ret;


}


//' @title testing copula GGM
//' @name copula
//' @export

// mixed data sampler
// [[Rcpp::export]]
Rcpp::List  copula(arma::mat z0_start,
                   arma::mat levels,
                   arma::vec K,
                   arma::mat Sigma_start,
                   int iter,
                   float delta,
                   float epsilon,
                   arma::vec idx) {

  // adapted from hoff 2008 for Bayesian hypothesis testing
  // with the matrix-F prior distribution for Theta

  // z0: latent data
  // levels: data matrix as sorted levels
  // K: levels in each columns

  // number of rows
  float n = z0_start.n_rows;

  // number of columns
  int k = z0_start.n_cols;

  // k by k identity mat
  arma::mat  I_k(k, k, arma::fill::eye);

  int nu = 1/ epsilon;
  // // #nu in Mulder & Pericchi (2018) formula (30) line 1.
  int nuMP = delta + k - 1 ;
  //
  // // #delta in Mulder & Pericchi (2018) formula (30) line 1.
  int deltaMP = nu - k + 1 ;

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
  arma::cube Theta_mcmc(k, k, iter, arma::fill::zeros);

  // partial correlations
  arma::mat pcors(k,k);
  arma::cube pcors_mcmc(k, k, iter, arma::fill::zeros);

  // correlations
  arma::mat  cors(k,k);
  arma::cube cors_mcmc(k, k, iter, arma::fill::zeros);

  // covariance matrix
  Sigma.slice(0) = Sigma_start;
  arma::cube Sigma_mcmc(k, k, iter, arma::fill::zeros);

  arma::vec lb(1);
  arma::vec ub(1);

  arma::mat mm(n,1);
  arma::mat ss(1,1);

  for(int  s = 1; s < iter; ++s){

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

    // correlation
    cors =  diagmat(1 / sqrt(Sigma.slice(0).diag())) *
      Sigma.slice(0) *
      diagmat(1 / sqrt(Sigma.slice(0).diag()));

    // partial correlations
    pcors = diagmat(1 / sqrt(Theta.slice(0).diag())) *
      Theta.slice(0) *
      diagmat(1 / sqrt(Theta.slice(0).diag()));

    pcors_mcmc.slice(s) =  -(pcors - I_k);
    cors_mcmc.slice(s) =  cors;
    Sigma_mcmc.slice(s) = Sigma.slice(0);
    Theta_mcmc.slice(s) = Theta.slice(0);
  }


  Rcpp::List ret;
  ret["pcors"] = pcors_mcmc;
  ret["cors"] =  cors_mcmc;
  ret["Theta"] = Theta_mcmc;
  ret["Sigma"] = Sigma_mcmc;
  return ret;
}




