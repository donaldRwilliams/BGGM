#' GGM comparison based on the predictive distribution
#'
#' @param Y_g1
#' @param Y_g2
#' @param Y_g3
#' @param Y_g4
#' @param contrast
#' @param groups
#'
#' @param samples numer of posterior and predictive samples
#' @return
#' @export
#'
#' @examples
GGM_compare_ppc <- function(Y_g1,
                            Y_g2,
                            Y_g3 = NULL,
                            Y_g4 = NULL,
                            contrast,
                            samples){
  # number of columns
  p <- ncol(Y_g1)

  # remove all but numbers from contrast
  contrast <- unlist(regmatches(contrast, gregexpr("[[:digit:]]+", contrast)))

  # ensure contrast is of length 2
  if(length(contrast) != 2) stop("contrast must be between two groups")

  # check for two groups
  if(is.null(Y_g3) && is.null(Y_g4)){
    # number of groups
    groups <- 2

    if(sum(as.numeric(contrast)) != 3){

      # check contrast is between groups 1 and 2
      stop("with two groups contrast must be between 1 and 2--e.g., 1 vs 2, 1 & 2, et.")

    }
  } else if(is.null(Y_g4) & !is.null(Y_g1) & !is.null(Y_g2)){
    # number of groups
    groups <- 3

    # check contrast
    check_contrast <-  any(colSums(c(contrast) == combn(1:groups, 2)) == 2)

    # if false notify of invalid contrast
    if(isFALSE(check_contrast)) stop("invalid contrast. please respecify")

  } else if(!is.null(Y_g4) & !is.null(Y_g3) & !is.null(Y_g2) & !is.null(Y_g1)){
    # number of groups
    groups <- 4

    # number of groups
    check_contrast <-  any(colSums(c(contrast) == combn(1:groups, 2)) == 2)

    # if false notify of invalid contrast
    if(isFALSE(check_contrast)) stop("invalid contrast. please respecify")

  } else {

    # invalid input (nothing recognized)
    stop("invalid input")

  }

  # contrast
  contrast <- paste("Y", contrast, sep = "_g", collapse = "_vs_")

  # 2 groups
  if(groups == 2){
    # group 1
    Y_g1 <- as.matrix(Y_g1)
    n1 <- nrow(Y_g1)

    # group 2
    Y_g2 <- as.matrix(Y_g2)
    n2 <- nrow(Y_g2)

    # scale the combined data
    Y_G <- scale(rbind(Y_g1, Y_g2), scale = T)

    # inverse scatter matrix
    S_G <- solve(t(Y_G) %*% Y_G)

    # M_0 posterior (group equality)
    post <- rWishart(samples, n1 + n2 - 1, S_G)


    if(contrast == "Y_g1_vs_Y_g2"){

      # predictive risk of the null model
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g2)) + 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g1))
    }
  }

  # 3 groups
  if(groups == 3){
    # group 1
    Y_g1 <- as.matrix(Y_g1)
    n1 <- nrow(Y_g1)

    # group 2
    Y_g2 <- as.matrix(Y_g2)
    n2 <- nrow(Y_g2)

    # group 3
    Y_g3 <- as.matrix(Y_g3)
    n3 <- nrow(Y_g3)

    # scale the combined data
    Y_G <- scale(rbind(Y_g1, Y_g2, Y_g3), scale = F)

    # inverse scatter matrix
    S_G <- solve(t(Y_G) %*% Y_G)

    # M_0 posterior (group equality)
    post <- rWishart(samples, n1 + n2 + n3 - 1, S_G)


    if(contrast == "Y_g1_vs_Y_g2"){

      # predictive risk of the null model
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g2)) + 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g1))
    }

    if(contrast == "Y_g1_vs_Y_g3"){

      # predictive risk of the null model
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n3, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g1))

    }

    if(contrast == "Y_g2_vs_Y_g3"){

      # predictive risk of the null model
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n2, n3, p))

      # observed risk
      obs_jsd <- 0.5 *  KL(unbiased_cov(Y_g2), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g2))

    }
  }

  # 4 groups
  if(groups == 4){
    # group 1
    Y_g1 <- as.matrix(Y_g1)
    n1 <- nrow(Y_g1)

    # group 2
    Y_g2 <- as.matrix(Y_g2)
    n2 <- nrow(Y_g2)

    # group 3
    Y_g3 <- as.matrix(Y_g3)
    n3 <- nrow(Y_g3)

    # group 4
    Y_g4 <- as.matrix(Y_g4)
    n4 <- nrow(Y_g4)

    # scale the combined data
    Y_G <- scale(rbind(Y_g1, Y_g2, Y_g3, Y_g4), scale = F)

    # inverse scatter matrix
    S_G <- solve(t(Y_G) %*% Y_G)

    # M_0 posterior (group equality)
    post <- rWishart(samples, n1 + n2 + n3 + n4 - 1, S_G)

    if(contrast == "Y_g1_vs_Y_g2"){

      # M_0 posterior (group equality)
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g2)) + 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g1))

    }

    if(contrast == "Y_g1_vs_Y_g3"){

      # M_0 posterior (group equality)
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n3, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g1))

    }

    if(contrast == "Y_g1_vs_Y_g4"){

      # M_0 posterior (group equality)
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n4, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g4)) + 0.5 * KL(unbiased_cov(Y_g4), unbiased_cov(Y_g1))

    }

    if(contrast == "Y_g2_vs_Y_g3"){

      # M_0 posterior (group equality)
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n2, n3, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g2))

    }

    if(contrast == "Y_g2_vs_Y_g4"){

      # M_0 posterior (group equality)
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n2, n4, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g4)) + 0.5 * KL(unbiased_cov(Y_g4), unbiased_cov(Y_g2))

    }

    if(contrast == "Y_g3_vs_Y_g4"){

      # M_0 posterior (group equality)
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n3, n4, p))

      # observed risk
      obs_jsd <- 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g4)) + 0.5 * KL(unbiased_cov(Y_g4), unbiased_cov(Y_g3))

    }
  }

  # compute p-value
  p_value_jsd <- mean(unlist(Mo_risk) > obs_jsd)

  # make data frame with loss and p-value
  p_value_df <- data.frame(loss = c("jsd"), p_value = c(p_value_jsd))

  # returned object
  returned <- list(p_value_df = p_value_df,
                   obs_jsd = obs_jsd,
                   contrast = contrast,
                   jsd_scores = unlist(Mo_risk))

  # assign class
  class(returned) <- "GGM_compare_ppc"

  returned
}

