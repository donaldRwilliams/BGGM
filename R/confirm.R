#' Confirmatory Hypothesis Testing
#'
#' @description Traditionally, Gaussian graphical models (GGM) are inherently exploratory. That is, automated model selection is performed. A key aspect of \strong{BGGM} is the ability to extend inference beyond exploratory and to
#' confirmatory hypothesis testing. This is accomplished by testing equality and/or inequality constraints for sets of edges (partial correlations).
#'
#' @param Y data matrix  (\emph{n} by  \emph{p}).
#' @param hypothesis hypothesis (or hypotheses) to be tested
#' @param prior_sd hypothesized standard deviation of the prior distribution
#' @param iter posterior and prior samples. 25,000 is the default, as it results in a more stable Bayes factor than using, say, 5,000.
#' @param cores number of cores for parallel computing. The default is 2, but this can be adjusted
#'
#' @return list of class \code{confirm}:
#'
#' \itemize{
#' \item \code{BF_matrix} matrix of Bayes factors for each hypothesis. Also includes the compliment
#' \item \code{post_prob} posterior hypothesis probabilities
#' \item \code{hypotheses} \code{hypothesis}
#' \item \code{call} match.call()
#' \item \code{p} number of variables
#' \item \code{n} number of observations
#' \item \code{iter} number of posterior samples
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to prior_sd)
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{returned_mats} contrast matrices
#' }
#'
#' @importFrom MASS ginv
#' @importFrom stats cov rbeta
#' @export
#'
#' @note Currently inequality and equality restrictions can be tested. The former is an ordering the respective edge sizes,
#' whereas the latter allows for testing whether certain edges are exactly the same.
#'
#' see \code{methods(class = "confirm")}
#'
#' @examples
#' \donttest{
#'
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # hypothesis
#' hypothesis <- c("1--2 > 1--3 > 1--4 > 1--5")
#'
#' # test inequality contraint
#' test_order <-  confirm(Y = Y, hypothesis  = hypothesis,
#'                       prior_sd = 0.5, iter = 50000,
#'                       cores = 2)
#' # summary
#' summary(test_order)
#'
#'
#'# test hypothesized directions
#'
#'# hypothesis
#'hypothesis <- c("(1--2, 1--3, 1--4)  <  0 < (1--6)")
#'
#'# test directions
#' test_directions <-  confirm(Y = Y, hypothesis  = hypothesis,
#'                       prior_sd = 0.5, iter = 50000,
#'                       cores = 2)
#'# summary
#' test_directions
#'}
confirm <- function(Y, hypothesis, prior_sd = 0.25,
                    iter = 25000,  cores = 2){

  # code taken and adapted from (with permission):
  # https://github.com/Jaeoc/lmhyp

  # set prior prob to 1--i.e., equal
  priorprob = 1

  returned_mats <-list()

  names_check <- unlist(strsplit( strsplit(hypothesis, " ")[[1]], "--"))

  names_check <- paste(names_check, collapse = " ")

  names_check <- unique(strsplit(gsub("[^[:alnum:] ]", "", names_check), " +")[[1]])

  if(any(names_check == "0")){
    names_check <- names_check[-which(names_check == "0" )]
  }




  if(any(names_check %in% colnames(Y))){
     if(any(grepl('[^[:alnum:]]', colnames(Y)))){
      stop("special characters not allowed in column names")
     }
    if(!all(names_check %in% colnames(Y))){
      stop("node names not found in the data")
    }

    hyp_temp <- BGGM:::convert_colnames(hyp = hypothesis, Y = Y)
  } else {

    hyp_temp <- hypothesis
  }

  # convert hypotheses
  hyp <- BGGM:::hyp_converter(hyp_temp)$hyp_converted

  # fit model
  fit <- explore(Y = Y, prior_sd = prior_sd,
                         iter = iter, cores = cores)

  # number of edges
  edges <- 0.5* fit$p * (fit$p -1)

  # posterior
  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:fit$cores, function(z) fit$samples[[z]]$fisher_z_post))[1:edges]

  # prior samples
  prior_samples <-     do.call(rbind.data.frame,
                               lapply(1:fit$cores, function(z)  fit$samples[[z]]$fisher_z_prior))[1:edges]

  for(loop in seq_along(hyp)){
    # remove spaces
    hyp2 <- gsub("[ \n]", "", hyp[[loop]])

    # check hypotheses
    if(!grepl("^[0-9a-zA-Z><=,;().-]+$", hyp2)) stop("Impermissable characters in hypotheses.")

    if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

    # split hypotheses
    hyp_out <- if(length(hyp) > 1) c("X < 0", "X = 0", "X > 0") else unlist(strsplit(hyp2, split = ";"))

    hypothesis <- gsub("[[:space:]]", "",
                       sub(pattern = "\n", "",
                           unlist(strsplit(hypothesis, split = ";"))))

    for(no in seq_along(hypothesis)){
      names(hypothesis)[no] <- paste0("H", no)
      names(hypothesis)[no] <- paste0("H", no)
    }


    # if(!isTRUE(priorprob == 1) && length(priorprob) < length(hyp_out)) stop("Use default equal priorprob or specify priorprob for all hypotheses")

    hyps <- unlist(strsplit(hyp2, split = ";"))

    BFu <- out_c_E <- out_f_E <- out_c_i_e <- out_f_i_e <- out_ci_lb <- out_ci_ub <- rep(NA, length = length(hyps))

    BFip_posterior <- if(any(!grepl("=", hyps))) {rep(NA, sum(!grepl("=", hyps)))} else{NULL}

    if(!is.null(BFip_posterior)) {

      R_i_all <- r_i_all <- BFi_draws_post <- BFi_draws_prior <- vector("list", length =  length(BFip_posterior))

      ineq_marker <- 0

      BFip_prior <- BFip_posterior
    }

    for(h in seq_along(hyps)){

      hyp2 <- hyps[[h]]

      # seperate
      step2 <- unlist(strsplit(hyp2, split = "[<>=,()]"))

      # which vars
      hyp_vars <- step2[grep("[a-zA-Z]+", step2)]



      if(sum(hyp_vars %in% colnames(posterior_samples)) != length(hyp_vars)) {
        stop("error in edge specification\n
             Hints:\n
             1) check that the first is smaller than the second (denoting the upper-triangular elements)-e.g., 1--2 and not 2--1\n
             2) alternatively, the variable might not exist--e.g., p = 10 but 1--11 (11 > 10)")
      }



      # check dupilicates
      if(any(duplicated(hyp_vars))) stop("Variables should occur only once in a hypothesis. Check semicolons.")

      framed <- BGGM:::framer(hyp2)

      mats <- BGGM:::create_matrices(framed = framed, varnames = colnames(posterior_samples))

      returned_mats[[h]] <- mats

      R_ei <- rbind(mats$R_e,mats$R_i)

      r_ei <- rbind(mats$r_e,mats$r_i)

      Rr_ei <- cbind(mats$R_ei,mats$r_ei)

      beta_zero <- MASS::ginv(mats$R_ei)%*%mats$r_ei



      if(mats$comparisons == "only equality"){

        Sigma0 <- cov(prior_samples)
        Sigma1 <- cov(posterior_samples)

        mu0 <- mats$R_e %*% colMeans(prior_samples)
        mu1 <- mats$R_e %*% colMeans(posterior_samples)


        s0 <- mats$R_e %*% Sigma0 %*% t(mats$R_e)
        s1 <- mats$R_e %*% Sigma1 %*% t(mats$R_e)
        #Hypothesis test
        log_BF <- mvtnorm::dmvnorm(x = t(mats$r_e), mean =   t(mu1), sigma = s1, log = TRUE) -
          mvtnorm::dmvnorm(x = t(mats$r_e), mean = t(mu0), sigma = s0, log = TRUE)

        f_E <- mvtnorm::dmvnorm(x = t(mats$r_e), mean =   t(mu1), sigma = s1, log = FALSE)
        c_E <- mvtnorm::dmvnorm(x = t(mats$r_e), mean = t(mu0), sigma = s0, log = TRUE)
        c_i_e <- f_i_e <- NA

        BF <- exp(log_BF)
        ci_lb <- ci_ub <- NA

      }

      if(mats$comparisons == "only inequality"){

        ineq_marker <- ineq_marker + 1

         if(Matrix::rankMatrix(mats$R_i)[[1]] == nrow(mats$R_i)){

          Sigma0 <- cov(prior_samples)
          Sigma1 <- cov(posterior_samples)

          mu0 <- mats$R_i %*% colMeans(prior_samples)
          mu1 <- mats$R_i %*% colMeans(posterior_samples)

          s0 <- mats$R_i %*% Sigma0 %*% t(mats$R_i)
          s1 <- mats$R_i %*% Sigma1 %*% t(mats$R_i)

          prior_prob <- mvtnorm::pmvnorm(lower = as.vector(mats$r_i),
                                         upper = Inf,
                                         mean = as.vector(mu0),
                                         sigma = s0)[1]

          posterior_prob <-  mvtnorm::pmvnorm(lower = as.vector(mats$r_i),
                                              upper = Inf,
                                              mean = as.vector(mu1),
                                              sigma = s1)[1]

          BF <- posterior_prob / prior_prob

          ci_lb <- ci_ub <- NA

          ci_draws_post <- ci_draws_prior <- NULL

         } else {

           mcrep <- 1e+06
           Sigma0 <- cov(prior_samples)
           Sigma1 <- cov(posterior_samples)

           draws_post <- mvtnorm::rmvt(n = mcrep,
                                       delta = colMeans(posterior_samples),
                                       sigma =  Sigma1)
           draws_prior <- mvtnorm::rmvt(n = mcrep, delta = mats$beta_zero,
                                        sigma = Sigma1)

           prior_prob <- mean(apply(draws_prior%*%t(mats$R_i) > rep(1, mcrep)%*%t(mats$r_i), 1, prod))
           posterior_prob <- mean(apply(draws_post%*%t(mats$R_i) > rep(1, mcrep)%*%t(mats$r_i), 1, prod))

           BF <- posterior_prob / prior_prob

             #Credibility interval
             x_post <- sum(apply(draws_post%*%t(mats$R_i) > rep(1, mcrep)%*%t(mats$r_i), 1, prod))
             x_prior <- sum(apply(draws_prior%*%t(mats$R_i) > rep(1, mcrep)%*%t(mats$r_i), 1, prod))
             ci_draws_post <- rbeta(1e4, x_post, 1 + mcrep - x_post)
             ci_draws_prior <- rbeta(1e4, x_prior, 1 + mcrep - x_prior)
             BF_draws <- ci_draws_post / ci_draws_prior

             ci_lb <- quantile(sort(BF_draws), 0.05)
             ci_ub <- quantile(sort(BF_draws), 0.95)

        }

        if(mats$comparisons == "both comparisons") {
          stop("inequality and equality constraints not currently supported")

        }



        c_i_e <- prior_prob
        f_i_e <- posterior_prob
        f_E <- c_E <- NA

        BFip_prior[ineq_marker] <- prior_prob
        BFip_posterior[ineq_marker] <- posterior_prob
        R_i_all[[ineq_marker]] <- mats$R_i
        r_i_all[[ineq_marker]] <- matrix(mats$r_i)
        BFi_draws_post[ineq_marker] <- list(ci_draws_post)
        BFi_draws_prior[ineq_marker] <- list(ci_draws_prior)
      }

      out_c_E[h] <- c_E
      out_f_E[h] <- f_E
      out_c_i_e[h] <- c_i_e
      out_f_i_e[h] <- f_i_e
      out_ci_lb[h] <- ci_lb
      out_ci_ub[h] <- ci_ub

      BFu[h] <- BF #hypothesis vs. unconstrained
      names(BFu)[[h]] <- paste0("H", h)

    }


  if(!is.null(BFip_posterior)){
      if(length(BFip_posterior) == 1){
        BFc <- (1 - BFip_posterior) / (1 - BFip_prior)
        comp_fie <- (1 - BFip_posterior)
        comp_cie <- (1 - BFip_prior)
        comp_fe <- comp_ce <- NA

        #Credibility interval
        if(length(Filter(Negate(is.null), BFi_draws_post)) == 0){
          BFc_ci_lb <- BFc_ci_ub <- NA
        }else{
          BFc_draws_post <- 1 - unlist(BFi_draws_post)
          BFc_draws_prior <- 1 - unlist(BFi_draws_prior)

          BFc_draws <- BFc_draws_post / BFc_draws_prior
          BFc_ci_lb <- quantile(sort(BFc_draws), 0.05)
          BFc_ci_ub <- quantile(sort(BFc_draws), 0.95)
        }
      } else{

        R_i_overlap <- do.call(rbind, R_i_all)
        r_i_overlap <- do.call(rbind, r_i_all)
        ineq_draws_prior <- mvnfast::rmvn(n = 1e4, mu = colMeans(prior_samples), sigma = cov(prior_samples))

        exhaustive <- mean(rowSums(ineq_draws_prior%*%t(R_i_overlap) > rep(1, 1e4)%*%t(r_i_overlap)) > 0)

        if(exhaustive == 1){
          BFc <- comp_cie <- comp_fie <- comp_fe <- comp_ce <- BFc_ci_lb <- BFc_ci_ub <- NULL
        } else{
          overlap <- mean(apply(ineq_draws_prior%*%t(R_i_overlap) > rep(1, 1e4)%*%t(r_i_overlap), 1, prod))

          if(overlap == 0){
            BFc <-  (1 - sum(BFip_posterior)) / (1 - sum(BFip_prior))
            comp_fie <- (1 - sum(BFip_posterior))
            comp_cie <- (1 - sum(BFip_prior))

            #Credibility interval
            if(length(Filter(Negate(is.null), BFi_draws_post)) == 0){
              BFc_ci_lb <- BFc_ci_ub <- NA
            }else{
              analytic_BFip_posterior <- BFip_posterior[which(unlist(lapply(BFi_draws_post, is.null)) == 1)]
              num_draws_BFip_posterior <- Filter(Negate(is.null),BFi_draws_post)
              BFc_draws_post <- 1 - sum(analytic_BFip_posterior) - Reduce(`-`, num_draws_BFip_posterior)

              analytic_BFip_prior <- BFip_prior[which(unlist(lapply(BFi_draws_prior, is.null)) == 1)]
              num_draws_BFip_prior <- Filter(Negate(is.null),BFi_draws_prior)
              BFc_draws_prior <- 1 - sum(analytic_BFip_prior) - Reduce(`-`, num_draws_BFip_prior)

              BFc_draws <- BFc_draws_post / BFc_draws_prior
              BFc_ci_lb <- quantile(sort(BFc_draws), 0.05)
              BFc_ci_ub <- quantile(sort(BFc_draws), 0.95)
            }
          } else{

            # ineq_draws_posterior <- mvtnorm::rmvt(n = 1e4, delta = betahat, sigma = vcov(object), df = n - k)

            ineq_draws_posterior <- mvnfast::rmvn(n = 1e4, colMeans(posterior_samples), sigma = cov(posterior_samples))

            constraints_prior <- Map(function(Ri, ri){apply(ineq_draws_prior%*%t(Ri) > rep(1,1e4)%*%t(ri), 1, prod)}, R_i_all, r_i_all)
            constraints_posterior <- Map(function(Ri, ri){apply(ineq_draws_posterior%*%t(Ri) > rep(1,1e4)%*%t(ri), 1, prod)}, R_i_all, r_i_all)

            prob_prior <- mean(Reduce(`+`, constraints_prior) > 0)
            prob_posterior <- mean(Reduce(`+`, constraints_posterior) > 0)

            BFc <- (1 - prob_posterior) / (1 - prob_prior)
            comp_fie <- (1 - prob_posterior)
            comp_cie <- (1 - prob_prior)

            #Credibility interval
            x_post <- 1e4 - sum(Reduce(`+`, constraints_posterior) > 0)
            x_prior <- 1e4 - sum(Reduce(`+`, constraints_prior) > 0)
            ci_draws_post <- rbeta(1e4, x_post, 1 + 1e4 - x_post)
            ci_draws_prior <- rbeta(1e4, x_prior, 1 + 1e4 - x_prior)
            BFc_draws <- ci_draws_post / ci_draws_prior

            BFc_ci_lb <- quantile(sort(BFc_draws), 0.05)
            BFc_ci_ub <- quantile(sort(BFc_draws), 0.95)
          }
          comp_fe <- comp_ce <- NA
        }
      }
    } else{
      BFc <- 1
      comp_ce <- 1
      comp_fe <- 1
      comp_fie <- comp_cie <- BFc_ci_lb <- BFc_ci_ub <- NA
    }

    if(!is.null(BFc)){names(BFc) <- "Hc"}
    BFu <- c(BFu, BFc)
    if(!isTRUE(priorprob == 1) && !is.null(BFc) && length(BFu) > length(priorprob)) {
      stop("Hypotheses have complement, add priorprob for it or use equal priorprobs.")
    }

    if(length(priorprob) > length(BFu)) stop("Too many priorprobs specified")
    if(!sum(priorprob) == 1){
      priorprob <- priorprob / sum(priorprob)
      warning(paste(c("priorprobs did not sum to 1 and have been normalized. Used priorprobs: ",
                      round(priorprob, 4)), collapse = " "))
    }
    out_hyp_prob <- BFu*priorprob / sum(BFu*priorprob)

    BF_matrix <- matrix(rep(BFu, length(BFu)), ncol = length(BFu), byrow = TRUE)
    BF_matrix <- t(BF_matrix / BFu)
    colnames(BF_matrix) <- rownames(BF_matrix) <- names(BFu)
    BF_matrix[is.nan(BF_matrix)] <- 0
    diag(BF_matrix) <- 1

    if(length(hyp) > 1){
      BF_matrices[[loop]] <- round(BF_matrix, digits = 3)
      list_post_prob[[loop]] <- out_hyp_prob
    }

  } #end exploratory loop


  if(length(hyp) > 1){
    matrix_post_prob <- matrix(do.call(rbind, list_post_prob), ncol = 3)

    colnames(matrix_post_prob) <- c("H1", "H2", "H3")
    varnames <- gsub("Intercept", "(Intercept)", varnames)
    names(BF_matrices) <- rownames(matrix_post_prob) <- varnames

    out <- list(BF_matrix = BF_matrices, post_prob = matrix_post_prob,
                hypotheses = hypothesis, BF_computation = "Not available when option 'exploratory' chosen.",
                BFu_CI = "Not available when option 'exploratory' chosen.")

  } else{
    out_c_i_e <- c(out_c_i_e, comp_cie)
    out_f_i_e <- c(out_f_i_e, comp_fie)
    out_c_E <- c(out_c_E, comp_ce)
    out_f_E <- c(out_f_E, comp_fe)
    out_ci_lb <- c(out_ci_lb, BFc_ci_lb)
    out_ci_ub <- c(out_ci_ub, BFc_ci_ub)

    BF_computation <- as.matrix(data.frame(out_c_E, out_c_i_e, c = out_c_E*out_c_i_e, out_f_E, out_f_i_e,
                                           f = out_f_E*out_f_i_e, BFu, PP_t = out_hyp_prob))
    colnames(BF_computation) <- c("c(E)", "C(I|E)", "c", "f(E)", "f(I|E)", "f", "B(t,u)", "PP(t)")

    BFu_ci <- as.matrix(data.frame(BFu, out_ci_lb, out_ci_ub))
    colnames(BFu_ci) <- c("B(t,u)", "lb. (5%)", "ub. (95%)")


    returned_object <- list(BF_matrix = round(BF_matrix, digits = 3),
                post_prob = out_hyp_prob,
                hypotheses = hypothesis,
                call = match.call(),
                p = fit$p, n = nrow(fit$dat),
                iter = fit$iter,
                delta = fit$delta,
                parcors_mat = fit$parcors_mat,
                returned_mats = returned_mats,
                BF_computation = BF_computation)
  }

  class(returned_object) <- c("BGGM", "confirm")
  returned_object

}
