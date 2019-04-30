#' GGM_compare_bf
#'
#' @param X1
#' @param X2
#' @param hyp
#' @param priorprob
#' @param cores
#' @param delta
#' @param nu
#' @param n_samples
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_bf <- function(Y_g1,
                           Y_g2,
                           Y_g3 = NULL,
                           Y_g4 = NULL,
                           hyp,
                           priorprob = 1,
                           cores,
                           delta,
                           nu,
                           n_samples){

  # check group number
  if(is.null(Y_g3) && is.null(Y_g4)){
    # number of groups
    groups <- 2

  } else if(is.null(Y_g4) & !is.null(Y_g1) & !is.null(Y_g2)){
    # number of groups
    groups <- 3

  } else if(!is.null(Y_g4) & !is.null(Y_g3) & !is.null(Y_g2) & !is.null(Y_g1)){
    # number of groups
    groups <- 4

  } else {
    # invalid input (nothing recognized)
    stop("invalid input")

  }

  p <- ncol(Y_g1)

  off_diag <-  0.5 * (p * (p - 1))


  # additions
  # Y_g1 <- MASS::mvrnorm(n = 300, rep(0, 16), Sigma = BGGM::ptsd_cor1, empirical = T)
  # Y_g2 <- MASS::mvrnorm(n = 950, rep(0, 16), Sigma = BGGM::ptsd_cor2, empirical = T)
  # Y_g3 <- MASS::mvrnorm(n = 950, rep(0, 16), Sigma = BGGM::ptsd_cor3, empirical = T)
  # nu = 20; delta = 20; n_samples = 5000; cores = 4
  #


  samples_1 <- BGGM:::sampling(as.matrix(Y_g1), nu = nu, delta = delta, n_samples = n_samples, cores)
  samples_2 <- BGGM:::sampling(as.matrix(Y_g2), nu = nu, delta = delta, n_samples = n_samples, cores)


  if(hyp == "all"){


    if(groups == 2) {


    post_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_post[,1:off_diag]))
    post_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_post[,1:off_diag]))


    prior_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_prior[,1:off_diag]))
    prior_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_prior[,1:off_diag]))


    post_diff <- post_1 - post_2
    prior_diff <- prior_1[,1] - prior_2[,1]

    post_dens <- apply(post_diff, 2, function(x){dnorm(0, mean(x), sd(x))})

    BF <-   post_dens  / dnorm(0, mean(prior_diff), sd(prior_diff))


    mat_temp <- BF_01 <- pcor_mean <- matrix(0, p, p)

    mat_temp[] <- unlist(lapply(1:p, function(x) paste(1:p, x, sep = "_")))

    BF_01[upper.tri(BF_01)] <- BF
    BF_01 <- BGGM:::symmteric_mat(BF_01)


    BF_10[upper.tri(BF_10)] <- 1 / BF
    BF_10 <- BGGM:::symmteric_mat(BF_10)


    prob_01 <- round(BF_01 / (1 + BF_01), 4)
    prob_10 <- round(1 - prob_01, 4)

    BF_01 <- round(BF_01, 4)
    BF_10 <- round(BF_10, 4)


    returned <- list(BF_01 = BF_01,
                     BF_10 = BF_10,
                     prob_01 = prob_01,
                     prob_10 = prob_10,
                     groups = groups)

    } # end groups == 2



    if(groups == 3){

    samples_3 <- BGGM:::sampling(as.matrix(Y_g3), nu = nu, delta = delta, n_samples = n_samples, cores)

    post_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_post[,1:off_diag]))
    post_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_post[,1:off_diag]))
    post_3 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_3[[x]]$fisher_z_post[,1:off_diag]))

    prior_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_prior[,1:off_diag]))
    prior_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_prior[,1:off_diag]))
    prior_3 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_3[[x]]$fisher_z_prior[,1:off_diag]))



    R <- matrix(c(1,-1,0,
                  0,1,-1),
                nrow = 2,
                byrow = T)


     for(i in 1:off_diag){

      mu_post <- c(mean(post_1[,i]), mean(post_2[,i]), mean(post_3[,i]) )
      cov_post <- cov(cbind(post_1[,i], post_2[,i], post_3[,i]))


      mu_prior <- c(mean(prior_1[,i]), mean(prior_2[,i]), mean(prior_3[,i]) )
      cov_prior <- cov(cbind(prior_1[,i], prior_2[,i], prior_3[,i]))


      mu1 <- R %*% mu_post
      s1 <- R %*% cov_post %*% t(R)


      mu0 <- rep(0, 2)
      s0 <- R %*% cov_prior %*% t(R)

      #Hypothesis test
      log_BF <- mvnfast::dmvn(X = rep(0,2), mu = mu1, sigma  = s1, log = TRUE) -
                mvnfast::dmvn(X = rep(0,2), mu = mu0, sigma = s0, log = TRUE)

       BF[i] <- exp(log_BF)
  }


    mat_temp <- BF_01 <- BF_10 <- matrix(0, p, p)
    mat_temp[] <- unlist(lapply(1:p, function(x) paste(1:p, x, sep = "_")))

    # contrast_name <- rep(c("1_vs_2", "1_vs_3", "2_vs_3"), each = off_diag)
    #
    # edge_name <- mat_temp[upper.tri(mat_temp)]
    # mean_diffs <- c(post_diff_1_vs_2, post_diff_1_vs_3, post_diff_2_vs_3)
    # sd_diffs <- c(sd_diff_1_vs_2, sd_diff_1_vs_3, sd_diff_2_vs_3)
    #
    # dat_1 <- cbind.data.frame(contrast = "1_vs_2",
    #                  edge_names,
    #                  mean_difference = post_diff_1_vs_2,
    #                  sd_difference = sd_diff_1_vs_2)
    #
    # dat_2 <- cbind.data.frame(contrast = "1_vs_3",
    #                           edge_names,
    #                           mean_difference = post_diff_1_vs_3,
    #                           sd_difference = sd_diff_1_vs_3)
    #
    # dat_3 <- cbind.data.frame(contrast = "2_vs_3",
    #                           edge_names,
    #                           mean_difference = post_diff_2_vs_3,
    #                           sd_difference = sd_diff_2_vs_3)
    #
    #
    # returned_list <- list(dat_1, dat_2, dat_3)
    # names(returned_list) <- c("1_vs_2", "1_vs_3", "2_vs_3")

   # BF <- round(BF, 4)

    BF_01[upper.tri(BF_01)] <- BF
    BF_01 <- BGGM:::symmteric_mat(BF_01)


    BF_10[upper.tri(BF_10)] <- 1 / BF
    BF_10 <- BGGM:::symmteric_mat(BF_10)


    prob_01 <- round(BF_01 / (1 + BF_01), 4)
    prob_10 <- round(1 - prob_01, 4)

    BF_01 <- round(BF_01, 4)
    BF_10 <- round(BF_10, 4)


    returned <- list(BF_01 = BF_01,
                     BF_10 = BF_10,
                     prob_01 = prob_01,
                     prob_10 = prob_10,
                     groups = groups)

   } # end groups three


    if(groups == 4){
      BF <- NA

      samples_3 <- BGGM:::sampling(as.matrix(Y_g3), nu = nu, delta = delta, n_samples = n_samples, cores)
      samples_4 <- BGGM:::sampling(as.matrix(Y_g4), nu = nu, delta = delta, n_samples = n_samples, cores)

      post_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_post[,1:off_diag]))
      post_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_post[,1:off_diag]))
      post_3 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_3[[x]]$fisher_z_post[,1:off_diag]))
      post_4 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_4[[x]]$fisher_z_post[,1:off_diag]))

      prior_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_prior[,1:off_diag]))
      prior_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_prior[,1:off_diag]))
      prior_3 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_3[[x]]$fisher_z_prior[,1:off_diag]))
      prior_4 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_4[[x]]$fisher_z_prior[,1:off_diag]))




      R <- matrix(c(1,-1,0, 0,
                    0,1,-1,0,
                    0,0,1,-1),
                  nrow = 3,
                  byrow = T)


      for(i in 1:off_diag){

        mu_post <- c(mean(post_1[,i]), mean(post_2[,i]), mean(post_3[,i]), mean(post_4[,i]))
        cov_post <- diag(diag(cov(cbind(post_1[,i], post_2[,i], post_3[,i], post_4[,i]))), 4)


        mu_prior <- c(mean(prior_1[,i]), mean(prior_2[,i]), mean(prior_3[,i]), mean(prior_4[,i]))
        cov_prior <- cov(cbind(prior_1[,i], prior_2[,i], prior_3[,i], prior_4[,i]))


        mu1 <- R %*% mu_post
        s1 <- R %*% cov_post %*% t(R)


        mu0 <- rep(0, 3)
        s0 <- R %*% cov_prior %*% t(R)

        #Hypothesis test
        log_BF <- mvnfast::dmvn(X = rep(0,3), mu = mu1, sigma  = s1, log = TRUE) -
                  mvnfast::dmvn(X = rep(0,3), mu = mu0, sigma = s0, log = TRUE)

        BF[i] <- exp(log_BF)
      }


      mat_temp <- BF_01 <- BF_10 <- matrix(0, p, p)
      # mat_temp[] <- unlist(lapply(1:p, function(x) paste(1:p, x, sep = "_")))
      #
      # contrast_name <- rep(c("1_vs_2",
      #                        "1_vs_3",
      #                        "1_vs_4",
      #                        "2_vs_3",
      #                        "2_vs_4",
      #                        "3_vs_4"), each = off_diag)
      #
      # edge_name <- mat_temp[upper.tri(mat_temp)]
      # mean_diffs <- c(post_diff_1_vs_2, post_diff_1_vs_3, post_diff_2_vs_3)
      # sd_diffs <- c(sd_diff_1_vs_2, sd_diff_1_vs_3, sd_diff_2_vs_3)
      #
      # dat_1 <- cbind.data.frame(contrast = "1_vs_2",
      #                           edge_names,
      #                           mean_difference = post_diff_1_vs_2,
      #                           sd_difference = sd_diff_1_vs_2)
      #
      # dat_2 <- cbind.data.frame(contrast = "1_vs_3",
      #                           edge_names,
      #                           mean_difference = post_diff_1_vs_3,
      #                           sd_difference = sd_diff_1_vs_3)
      #
      # dat_3 <- cbind.data.frame(contrast = "2_vs_3",
      #                           edge_names,
      #                           mean_difference = post_diff_2_vs_3,
      #                           sd_difference = sd_diff_2_vs_3)
      #
      #
      # returned_list <- list(dat_1, dat_2, dat_3)
      # names(returned_list) <- c("1_vs_2", "1_vs_3", "2_vs_3")

      # BF <- round(BF, 4)

      BF_01[upper.tri(BF_01)] <- BF
      BF_01 <- BGGM:::symmteric_mat(BF_01)


      BF_10[upper.tri(BF_10)] <- 1 / BF
      BF_10 <- BGGM:::symmteric_mat(BF_10)


      prob_01 <- round(BF_01 / (1 + BF_01), 4)
      prob_10 <- round(1 - prob_01, 4)

      BF_01 <- round(BF_01, 4)
      BF_10 <- round(BF_10, 4)


      returned <- list(BF_01 = BF_01,
                       BF_10 = BF_10,
                       prob_01 = prob_01,
                       prob_10 = prob_10,
                       groups = groups)

    }

   returned

}

if(hyp != "all"){

    hyp_new <- stringr::str_remove_all(hyp, "_")

    # gsub("(\\.+|[[:punct:]])", " \\1 ", hyp)

    pcor_names <- BGGM:::pcor_name_helper(hyp_new)

    gsub("[.][A-z]", ", [A-z]", hyp)

    pcor_names <- pcor_names[grep("[a-zA-Z]+", pcor_names)]



    post_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_post[,pcor_names]))
    post_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_post[,pcor_names]))


    prior_1 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_1[[x]]$fisher_z_prior[,pcor_names]))
    prior_2 <- do.call(rbind.data.frame, lapply(1:cores, function(x)    samples_2[[x]]$fisher_z_prior[,pcor_names]))


    post_diff <- post_1 - post_2
    prior_diff <- prior_1 - prior_2

    for(loop in seq_along(hyp_new)){

      # remove spaces
      hyp2 <- gsub("[ \n]", "", hyp_new[[loop]])
      #hyp2 <- gsub("(\\(Intercept\\))", "Intercept", hyp2)

      # check for permissable characters
      if(!grepl("^[0-9a-zA-Z><=,;().-]+$", hyp2)) stop("Impermissable characters in hypotheses.")

      # check to ensure there no combined comparison
      if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

      # make varnames individual characters
      step1 <- unlist(strsplit(hyp2, split = "[<>=,;()]"))

      # remove all but varnames
      input_vars <- step1[grep("[a-zA-Z]+", step1)]

      # check varnames are in the input vars
      # note: names are numbers spelled out in BGGM
      if(!all(input_vars %in% pcor_names)) stop("Hypothesis variable(s) not in object, check spelling")

      # split string at ";"
      hyp_out <- unlist(strsplit(hyp, split = ";"))

      # name the hypotheses H_{i}
      for(no in seq_along(hyp_out)){

        names(hyp_out)[no] <- paste0("H", no)

      }

      # check custom prior probabilities
      if(!isTRUE(priorprob == 1) && length(priorprob) < length(hyp_out)) {

        stop("Use default equal priorprob or specify priorprob for all hypotheses")

      }

      # split string at ";"
      hyps <- unlist(strsplit(hyp2, split = ";"))

      # storage
      BFu <-
        out_c_E <-
        out_f_E <-
        out_c_i_e <-
        out_f_i_e <-
        out_ci_lb <-
        out_ci_ub <- rep(NA, length = length(hyps))

      # number of (in)equality hypotheses
      BFip_posterior <- if(any(!grepl("=", hyps))) {

        rep(NA, sum(!grepl("=", hyps)))

      } else{ NULL

      }

      # if any inequality hypotheses
      if(!is.null(BFip_posterior)) {

        # storage for inequality
        R_i_all <-
          r_i_all <-
          BFi_draws_post <-
          BFi_draws_prior <- vector("list", length =  length(BFip_posterior))

        # marker for inequality
        ineq_marker <- 0

        # prior of same length
        BFip_prior <- BFip_posterior

      }



      ###########################
      #### end of 1st loop ######
      ###########################

      ############################################
      ################ second loop ###############
      ############################################
      # evaluate each hypotheses

      for(h in seq_along(hyps)){

        # temp hypothesis
        hyp2 <- hyps[[h]]

        # remove all but numbers and letters
        step2 <- unlist(strsplit(hyp2, split = "[<>=,()]"))

        # remove all all numbers
        hyp_vars <- step2[grep("[a-zA-Z]+", step2)]

        # check for duplicates
        if(any(duplicated(hyp_vars))) stop("Variables should occur only once in a hypothesis. Check semicolons.")


        # framing for matrices
        framed <- BGGM:::framer(hyp2)

        # matrix info
        matrix_info <- create_matrices(framed, varnames = pcor_names)

        comparisons <- matrix_info$comparisons


        # only equality hypotheses
        if(comparisons == "only equality"){

          # retrieve matrix info from function
          R_e <- matrix_info$R_e
          r_e <- matrix_info$r_e

          beta_zero <- matrix_info$beta_zero



          Sigma0 <- cov(prior_diff)
          Sigma1 <- cov(post_diff)

          mu0 <- R_e %*% beta_zero
          mu1 <- R_e %*% colMeans(post_diff)


          s0 <- R_e %*% Sigma0 %*% t(R_e)
          s1 <- R_e %*% Sigma1 %*% t(R_e)


          #Hypothesis test
          log_BF <- mvtnorm::dmvnorm(x = t(r_e), mean = mu1, sigma  = s1, log = TRUE) -
                    mvtnorm::dmvnorm(x = t(r_e), mean = mu0, sigma = s0, log = TRUE)

          f_E <- mvtnorm::dmvnorm(x = t(r_e), mean = mu1, sigma  = s1, log = FALSE)
          c_E <- mvtnorm::dmvnorm(x = t(r_e), mean = mu0, sigma = s0, log = FALSE)


          c_i_e <- f_i_e <- NA

          BF <- exp(log_BF)

          ci_lb <- ci_ub <- NA
        }


        out_c_E[h] <- c_E
        out_f_E[h] <- f_E
        out_c_i_e[h] <- c_i_e
        out_f_i_e[h] <- f_i_e
        out_ci_lb[h] <- ci_lb
        out_ci_ub[h] <- ci_ub

        BFu[h] <- BF #hypothesis vs. unconstrained
        names(BFu)[h] <- paste0("H", h)

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

          ineq_draws_prior <- MASS::mvrnorm(n = 1e4, mu =   matrix_info$beta_zero, Sigma = cov(prior))
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
              ineq_draws_posterior <- MASS::mvrnorm(n = 1e4, mu  = colMeans(post), Sigma = cov(post))

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


    }  #end exploratory loop

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

    returned <- list(BF_matrix = round(BF_matrix, digits = 3),
                     post_prob = out_hyp_prob,
                     hypotheses = hyp_out,
                     BF_computation = BF_computation,
                     BFu_CI = BFu_ci)

  }
  class(returned) <- "GGM_compare_bf"
  returned
}



