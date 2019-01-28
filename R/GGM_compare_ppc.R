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
GGM_compare_ppc <- function(Y_g1, Y_g2, Y_g3 = NULL, Y_g4 = NULL,
                            contrast, groups, samples) {

  # number of columns
  p <- ncol(Y_g1)



  if(groups == 2){
    Y_g1 <- as.matrix(Y_g1)
    Y_g2 <- as.matrix(Y_g2)


    n1 <- nrow(Y_g1)
    n2 <- nrow(Y_g2)

    Y_G <- scale(rbind(Y_g1, Y_g2), scale = F)
    S_G <- solve(t(Y_G) %*% Y_G)

    post <- rWishart(samples, n1 + n2 - 1, S_G)


    if(contrast == "Y_g1_vs_Y_g2"){

      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))

      obs_jd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g2)) + 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g1))

      pc1 <- cov2cor(unbiased_cov(Y_g1)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_g2)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      obs_sse <- sum((pc1 - pc2)^2)

    }

  }

  if(groups == 3){
    Y_g1 <- as.matrix(Y_g1)
    Y_g2 <- as.matrix(Y_g2)
    Y_g3 <- as.matrix(Y_g3)

    n1 <- nrow(Y_g1)
    n2 <- nrow(Y_g2)
    n3 <- nrow(Y_g3)

    Y_G <- scale(rbind(Y_g1, Y_g2, Y_g3), scale = F)
    S_G <- solve(t(Y_G) %*% Y_G)
    post <- rWishart(samples, n1 + n2 + n3 - 1, S_G)


    if(contrasts == "Y_g1_vs_Y_g2"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))

      obs_jd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g2)) + 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g1))

      pc1 <- cov2cor(unbiased_cov(Y_rep2)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep1)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)

    }

    if(contrasts == "Y_g1_vs_Y_g3"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n3, p))

      obs_jd <- 0.5 * BGGM::KL(mle(Y_g1), mle(Y_g3)) +
                0.5 * BGGM::KL(mle(Y_g3), mle(Y_g1))

      pc1 <- cov2cor(unbiased_cov(Y_rep3)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep1)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)


    }

   if(contrasts == "Y_g2_vs_Y_g3"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n2, n3, p))

      obs_jd <- 0.5 *  KL(unbiased_cov(Y_g2), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g2))

      pc1 <- cov2cor(unbiased_cov(Y_rep3)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep2)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)

    }

  }
  if(groups == 4){
    Y_g1 <- as.matrix(Y_g1)
    Y_g2 <- as.matrix(Y_g2)
    Y_g3 <- as.matrix(Y_g3)
    Y_g4 <- as.matrix(Y_g4)

    n1 <- nrow(Y_g1)
    n2 <- nrow(Y_g2)
    n3 <- nrow(Y_g3)
    n4 <- nrow(Y_g4)

    Y_G <- scale(rbind(Y_g1, Y_g2, Y_g3, Y_g4), scale = F)
    S_G <- solve(t(Y_G) %*% Y_G)

    post <- rWishart(samples, n1 + n2 + n3 + n4 - 1, S_G)

    if(contrasts == "Y_g1_vs_Y_g2"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))

      obs_jd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g2)) + 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g1))

      pc1 <- cov2cor(unbiased_cov(Y_rep2)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep1)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)


    }

    if(contrasts == "Y_g1_vs_Y_g3"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n3, p))


      obs_jd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g1))

      pc1 <- cov2cor(unbiased_cov(Y_rep3)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep1)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)



    }

    if(contrasts == "Y_g1_vs_Y_g4"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n1, n4, p))

      obs_jd <- 0.5 * KL(unbiased_cov(Y_g1), unbiased_cov(Y_g4)) + 0.5 * KL(unbiased_cov(Y_g4), unbiased_cov(Y_g1))

      pc1 <- cov2cor(unbiased_cov(Y_rep4)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep1)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)


    }

    if(contrasts == "Y_g2_vs_Y_g3"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n2, n3, p))


      obs_jd <- 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g3)) + 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g2))

      pc1 <- cov2cor(unbiased_cov(Y_rep3)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep2)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)


    }

    if(contrasts == "Y_g2_vs_Y_g4"){
      Mo_risk <-  lapply(1:samples, function(x) BGGM:::Mo_risk_help(x, post, n2, n4, p))

      obs_jd <- 0.5 * KL(unbiased_cov(Y_g2), unbiased_cov(Y_g4)) + 0.5 * KL(unbiased_cov(Y_g4), unbiased_cov(Y_g2))

      pc1 <- cov2cor(unbiased_cov(Y_rep4)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep2)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)

    }

    if(contrasts == "Y_g3_vs_Y_g4"){

      Mo_risk <-  lapply(1:samples, function(x) Mo_risk_help(x, post, n3, n4, p))

      obs_jd <- 0.5 * KL(unbiased_cov(Y_g3), unbiased_cov(Y_g4)) + 0.5 * KL(unbiased_cov(Y_g4), unbiased_cov(Y_g3))

      pc1 <- cov2cor(unbiased_cov(Y_rep4)) * -1
      pc1 <- pc1[upper.tri(pc1)]

      pc2 <- cov2cor(unbiased_cov(Y_rep3)) * -1
      pc2 <- pc2[upper.tri(pc2)] * -1

      sse <- sum((pc1 - pc2)^2)
    }


  }

  scores <- do.call(rbind.data.frame, Mo_risk)

  jd_scores <- subset(scores, loss == "jd")
  sse_scores <- subset(scores, loss == "sse")


  p_value_jd <- mean(jd_scores$score > obs_jd)
  p_value_sse <- mean(sse_scores$score > obs_sse)

  p_value_df <- data.frame(loss = c("jd", "sse"), p_value = c(p_value_jd, p_value_sse))

  returned <- list(p_value_df = p_value_df, scores = scores)
  returned


}

