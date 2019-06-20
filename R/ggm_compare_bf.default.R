#' Compare GGMs with Bayesian Hypothesis Testing
#'
#' @description Compare GGMs with the Bayes factor. This method allows for
#' assessing (relative) evidence for edge equality or edges differences across any number of groups. Further, confirmatory hypothesis testing
#' can be used to test predictions or expectations regarding difference or similiarities in different groups (e.g., male vs. female).
#'
#' @param ... data matrices (\emph{n} \times  \emph{p}). Requires at least two.
#' @param prior_sd hypothesized standard deviation for the edges or partial correlations.
#' @param iter number of posterior samples.
#' @param cores number of cores for parallel computing. The default is 2, but this can be adjusted.
#' @param hypothesis NULL
#'
#' @return list of class \code{ggm_compare_bf}
#' \itemize{
#' \item \code{BF_01} Bayes factors in favor of the null hypothesis
#' \item \code{p} number of variables
#' \item  \code{info} list of information about the data matrices
#'
#' \itemize{
#' \item \code{dat} list containing the data matrices
#' \item \code{dat_info} sample size for each data matrix
#' \item \code{pairwise} matrix of pairwise combinations
#' }
#' \item \code{iter} number of posterior and prior samples
#' \item \code{call} match.call()
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to \code{prior_sd})
#' }
#'
#' @note After fitting, use \code{select} to determine which partial correlations were different or the same (i.e., evidence for the null hypothesis of
#' equality). This assumes \code{hypothesis = NULL}. If a hypothesis is tested, then use \code{summary} which provides
#' information including the Bayes factors and posterior probabilities for each hypothesis.
#'
#' @export
#'
#'
#' @examples
#'
#'# assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' bf_ggm <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .5, iter = 2500, cores = 2)
#' select( bf_ggm, BF_cut = 3)

ggm_compare_bf.default <- function(..., prior_sd, hypothesis = NULL, iter = 25000, cores = 2){

  priorprob <- 1
  info <- BGGM:::Y_combine(...)

  groups <- length(info$dat)

  delta <- BGGM:::delta_solve(prior_sd)

  fits <- lapply(1:groups, function(x)   BGGM:::sampling(info$dat[[x]], nu = 10000,
                                                         delta = delta,
                                                         n_samples = iter,
                                                         cores = cores))

  # p
  p <- info$dat_info$p[1]

  edges <- 0.5 * (p * (p -1))

  prior_samps <- list()
  post_samps <- list()

  for(i in 1:groups){

    prior_samps[[i]] <- do.call(rbind.data.frame, lapply(1:cores, function(x)   fits[[i]][[x]]$fisher_z_prior[,1:edges]))

    post_samps[[i]] <- do.call(rbind.data.frame, lapply(1:cores, function(x)   fits[[i]][[x]]$fisher_z_post[,1:edges]))

  }


  if(is.null(hypothesis)){

    post_string <- list()
    prior_string <- list()

    for(i in 1:groups){

      post_string[[i]] <-   paste0("post_samps[[", i, "]][,",  1:edges, "]", sep = "")

      prior_string[[i]] <-   paste0("prior_samps[[", i, "]][,",  1:edges, "]", sep = "")

    }

    groups_as_words <- BGGM:::numbers2words(1:groups)

    hyp <- paste( groups_as_words, sep = " ", collapse = "=")

    framed <- BGGM:::framer(hyp)

    mats <- BGGM:::create_matrices(framed = framed, varnames = groups_as_words)


    post_string <- do.call(rbind, post_string)
    prior_string <- do.call(rbind, prior_string)

    BF <- NA
    for(i in 1:edges){

      # temporary string
      temp_string <- paste("cov(cbind(", paste(post_string[,i], sep = " ", collapse = ","), "))")

      # posterior covariance
      cov_post <- eval(parse( text = temp_string))

      # temporary string
      temp_string <- paste("colMeans(cbind(", paste(post_string[,i], sep = " ", collapse = ","), "))")

      # posterior mean
      post_mean <- eval(parse( text = temp_string))

      # temporary string
      temp_string <- paste("cov(cbind(", paste(prior_string[,i], sep = " ", collapse = ","), "))")

      # prior covariance
      cov_prior <- eval(parse( text = temp_string))

      # tranformed posterior
      mu1 <- mats$R_e %*% post_mean
      s1 <- mats$R_e %*% cov_post %*% t(mats$R_e)

      # transformed prior
      mu0 <- mats$R_e %*% rep(0, groups)
      s0 <- mats$R_e %*% cov_prior %*% t(mats$R_e)

      # bayes factor
      log_BF <- mvnfast::dmvn(X = t(mats$r_e), mu = mu1, sigma = s1, log = TRUE) -
        mvnfast::dmvn(X = t(mats$r_e), mu = mu0, sigma = s0, log = TRUE)

      BF[i] <- exp(log_BF)
    }
    BF_01 <- matrix(0, p, p)

    BF_01[upper.tri(BF_01)] <- BF

    BF_01 <- BGGM:::symmteric_mat(BF_01)

    returned_object <- list(BF_01 = BF_01,
                            p = p,
                            info = info,
                            iter = iter,
                            call = match.call(),
                            delta = delta)

  }

  if(!is.null(hypothesis)){


    temp_names <- lapply(1:groups, function(x) paste0("g", BGGM:::numbers2words(x),  colnames(prior_samps[[1]])))

    posterior_samples <- do.call(cbind.data.frame, post_samps)
    colnames(posterior_samples) <- unlist(temp_names)


    prior_samples <- do.call(cbind.data.frame, prior_samps)
    colnames(prior_samples) <- unlist(temp_names)

    returned_mats <-list()

    hyp <- BGGM:::hyp_converter(hypothesis)$hyp_converted
    hyp <- gsub( hyp, pattern = "_", replacement = "")

    for(loop in seq_along(hyp)){
      # remove spaces
      hyp2 <- gsub("[ \n]", "", hyp[[loop]])

      # check hypotheses
      if(!grepl("^[0-9a-zA-Z><=,;().-]+$", hyp2)) stop("Impermissable characters in hypotheses.")

      if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

      # split hypotheses
      hyp_out <- if(length(hyp) > 1) c("X < 0", "X = 0", "X > 0") else unlist(strsplit(hyp2, split = ";"))

      hypothesis <- gsub("[[:space:]]", "",  sub(pattern = "\n", "", unlist(strsplit(hypothesis, split = ";"))))

      for(no in seq_along(hyp_out)){
        names(hyp_out)[no] <- paste0("H", no)
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

        # # extract posterior samples
        # post_hyp <- posterior_samples[, hyp_vars]
        #
        # # extract prior samples
        # prior_hyp <- prior_samples[, hyp_vars]


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

          } else { # hyp mat not full rank

            Sigma0 <- cov(prior_samples)
            Sigma1 <- cov(posterior_samples)

            draws_post <- mvnfast::rmvn(1e+06, mu = colMeans(posterior_samples), sigma = Sigma1)
            draws_prior <- mvnfast::rmvn(1e+06, mu = mats$beta_zero, sigma = Sigma0)


            prior_prob <- mean(apply(draws_prior%*%t(mats$R_i) > rep(1, 1e+06)%*%t(mats$r_i), 1, prod))
            posterior_prob <- mean(apply(draws_post%*%t(mats$R_i) > rep(1, 1e+06)%*%t(mats$r_i), 1, prod))

            BF <- posterior_prob / prior_prob

            #Credibility interval
            x_post <- sum(apply(draws_post%*%t(mats$R_i) > rep(1, 1e+06)%*%t(mats$r_i), 1, prod))
            x_prior <- sum(apply(draws_prior%*%t(mats$R_i) > rep(1, 1e+06)%*%t(mats$r_i), 1, prod))
            ci_draws_post <- rbeta(1e4, x_post, 1 + 1e+06 - x_post)
            ci_draws_prior <- rbeta(1e4, x_prior, 1 + 1e+06 - x_prior)
            BF_draws <- ci_draws_post / ci_draws_prior

            ci_lb <- quantile(sort(BF_draws), 0.05)
            ci_ub <- quantile(sort(BF_draws), 0.95)
          } # end of else

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

          # ineq_draws_prior <- mvtnorm::rmvt(n = 1e4, delta = colMeans(prior_samples), sigma = vcov(object) * (n - k) / (n*b - k), df = (n*b - k))
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
                              BF_computation = BF_computation,
                              BFu_CI = BFu_ci,
                              call = match.call(),
                              p = p, info = info,
                              iter = iter,
                              delta = delta,
                              returned_mats = returned_mats)
    }
    }
  class(returned_object) <- "ggm_compare_bf"
  returned_object
}
