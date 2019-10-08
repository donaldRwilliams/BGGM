#' Select Graphical Structure with Estimation Based Methods
#' @name select.estimate
#' @description This allows for not only estimating the conditional dependence structure, that is non-zero edges, but also the conditional **in**dependence
#' structure (evidence for no relation). For the latter, the region of practical equivalence must be specified
#'
#' @param object object of class \code{estimate.default}
#' @param cred credible interval width used for the decision rule
#' @param rope region of practical equivalence
#' @param prob posterior probability (see notes)
#' @param ... not currently used
#'
#' @return An object of class \code{select.estimate}:
#'
#' \code{analytic = TRUE}:
#' \itemize{
#' \item \code{partials_non_zero} selected partial correlation matrix
#' \item \code{adjacency_non_zero} adjacency matrix for the selected edges
#' \item \code{ci} credible interval width
#' \item \code{analytic} TRUE
#' \item \code{pcors_samples} posterior samples
#' }
#'
#'\code{analytic = FALSE}:
#'\itemize{
#' \item \code{partials_non_zero} selected partial correlation matrix
#' \item \code{adjacency_non_zero} adjacency matrix for the selected edges
#' \item \code{pcor_sd} posterior standard deviation
#' \item \code{ci} credible interval width
#' \item \code{rope} NULL
#'
#'}
#'
#' \code{credible interval}:
#' \itemize{
#' \item \code{partials_non_zero} selected partial correlation matrix  (outside of the rope)
#' \item \code{adjacency_non_zero} adjacency matrix for the selected edges (outside of the rope)
#' \item \code{partials_zero} partials in the rope
#' \item \code{adjaceny_zero} adjacency in the rope
#' \item \code{pcor_sd} posterior standard deviation
#' \item \code{call} match.call()
#' \item \code{rope} specified rope
#' \item \code{in_rope} probability in the rope
#' \item \code{pcors_samples} posterior samples
#' }
#'
#'
#' @note The region of practical equivalence allows for assessing whether an edge is practically zero. In other words, conditional independence
#' (\eqn{\rho = 0}). The argument \code{prob} is then the posterior probability that must be in (practically zero edges) and out
#'  (practically zero edges) of the rope.
#'
#' @examples
#'
#' # Analytic = TRUE
#'# p = 5
#' Y <- BGGM::bfi[,1:5]
#'
#' # analytic solution
#' fit_analytic <- estimate(Y, analytic = TRUE)
#'
#' # select E
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' # non-zero partial correlations
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$adjacency_non_zero
#' @export
select.estimate <- function(object, cred = 0.95, rope = NULL, prob = 0.95, ...){

  x <- object
  ci_width <- cred

  # check object class
  if(class(x) !=  "estimate"){
    stop("Must be an object class bayes_estimate")
  }

  # ensure ci_width is allowed
  if(ci_width >= 1 | ci_width <= 0){
    stop("ci_width must be between 0 and 1")
  }

  if(isFALSE(x$analytic) ){

  pcor_samples <- x$posterior_samples[,  grep("pcors", colnames(x$posterior_samples))]
  pcor_sd <- matrix(0, x$p, x$p)
  pcor_sd[] <- apply(pcor_samples, 2, sd)

  if(is.null(rope)){
    # matrices for selected edges and the adjacency matrix
    adjacency_mat <- matrix(0, x$p, x$p)

    adjacency_mat[] <- apply(pcor_samples, 2, ci_helper, ci_width)

    returned_object <- list(partials_non_zero = x$parcors_mat * adjacency_mat,
                            pcor_sd = pcor_sd,
                            adjacency_non_zero = adjacency_mat,
                            call = match.call(),
                            ci = ci_width,
                            rope = rope,
                            pcor_samples  = pcor_samples)
  }

  if(is.numeric(rope)){
    message("cred is ignored")
    if(!is.numeric(prob)){
      stop("prob must be specificed (0 - 1) when rope = TRUE")
    }
    in_rope <- apply(pcor_samples, 2, rope_helper,  rope )

    out_rope <- 1 - in_rope
    nonzero_mat <- zero_mat <- matrix(0, x$p, x$p)
    nonzero_mat[] <- ifelse(out_rope >= prob, 1, 0)
    zero_mat[] <- ifelse(in_rope >= prob, 1, 0)

    returned_object <- list(partials_non_zero = x$parcors_mat * nonzero_mat,
                            adjacency_non_zero = nonzero_mat,
                            partial_zero = x$parcors_mat * zero_mat,
                            adjacency_zero = zero_mat,
                            pcor_sd = pcor_sd,
                            call = match.call(),
                            rope = rope,
                            prob = prob,
                            in_rope = in_rope,
                            pcor_samples = pcor_samples)
  }
  }
  if(!isFALSE(x$analytic)){
    crit <- abs(qnorm((1 - ci_width) /2))

    mat_sig <- matrix(0, x$p, x$p)
    up <-  x$fit$inv_mu[upper.tri(x$fit$inv_mu)] +   sqrt(x$fit$inv_var[upper.tri(x$fit$inv_var)]) * crit
    low <-  x$fit$inv_mu[upper.tri(x$fit$inv_mu)] -  sqrt(x$fit$inv_var[upper.tri(x$fit$inv_var)]) * crit

    mat_sig[upper.tri(mat_sig)] <- ifelse(low < 0 & up > 0, 0, 1)
    mat_sig <- symmteric_mat(mat_sig)

   returned_object <- list(partials_non_zero = x$fit$partial * mat_sig,
                           adjacency_non_zero = mat_sig,
                           call = match.call(),
                           ci = ci_width,
                           analytic = x$analytic)

  }
  class(returned_object) <- "select.estimate"
  returned_object$call <- match.call()
  returned_object
}

#' @title S3 select method
#' @name select
#' @description S3 select method
#' @param object object of class \code{estimate}, \code{explore}, or ..
#' @param ... not currently used
#' @return \code{select} works with the following methods:
#' \itemize{
#' \item \code{\link{select.estimate}}
#' \item \code{\link{select.explore}}
#' \item \code{\link{select.ggm_compare_estimate}}
#' }
#' @export
select <- function(object,...){
  UseMethod("select", object)
}

#' @name print.select.estimate
#' @title  Print method for \code{select.estimate} objects
#'
#' @param x An object of class \code{select.estimate}
#'
#' @param ... currently ignored
#' @seealso \code{\link{select.estimate}}
#'
#' @export
print.select.estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.numeric(x$rope)){
    cat("Type: Selected Graph (Sampling) \n")
  } else{
    cat("Type: Selected Graph (Analytic Solution) \n")
}
  if(is.null(x$rope)){
    cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
    cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    } else{
    cat("Probability:", x$prob, "\n")
    cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
    cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    }
  }

#' @name summary.select.estimate
#' @title Summary method for \code{select.estimate} objects
#'
#' @param object An object of class \code{select.estimate}
#' @param summarize if \code{TRUE} partial correlations and credible intervals are provided
#' @param ... currently ignored
#' @seealso \code{\link{select.estimate}}
#' @export
summary.select.estimate <- function(object, summarize = FALSE, ...){
  x <- object
  p <- ncol(x$partials_non_zero)
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(!is.null(x$analytic)){
    cat("Type: Selected Graph (Analytic Solution) \n")
  } else{
    cat("Type: Selected Graph (Sampling) \n")

  }

  if(isFALSE(summarize)){
    if(is.null(x$rope)){
      cat("Credible Interval:",  gsub("*0.","", formatC( round(x$ci, 4), format='f', digits=2)), "% \n")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Selected:\n \n")
      colnames( x$partials_non_zero)  <- 1:p
      row.names( x$partials_non_zero) <- 1:p
      colnames( x$partials_non_zero)  <- 1:p
      row.names( x$partials_non_zero) <- 1:p
      cat("Partial correlations \n \n")
      print(x$partials_non_zero, digits = 2)
      cat("--- \n \n")
      cat("Adjacency \n \n")
      colnames(x$adjacency_non_zero) <- 1:p
      row.names(x$adjacency_non_zero) <- 1:p
      print(x$adjacency_non_zero)
      cat("--- \n")

    } else{
      cat("Probability:", x$prob, "\n")
      cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Selected:\n \n")
      colnames(x$partials_non_zero) <- 1:p
      row.names(x$partials_non_zero) <- 1:p
      cat("Partial correlations \n \n")
      print(x$partials_non_zero, digits = 2)
      cat("--- \n \n")
      cat("Adjacency non-zero \n \n")
      colnames(x$adjacency_non_zero) <- 1:p
      rownames(x$adjacency_non_zero) <- 1:p
      print(x$adjacency_non_zero)
      cat("--- \n \n")
      cat("Adjacency zero \n \n")
      colnames(x$adjacency_zero) <- 1:p
      rownames(x$adjacency_zero) <- 1:p
      print(x$adjacency_zero)
    }
  }
  if(isTRUE(summarize)){
    if(isTRUE(x$analytic)){
      stop("summary not available for the analytic solution")
    }
    if(is.null(x$rope)){
      p <- ncol(x$partials_non_zero)
      mat_names <- mu_mat <- ci_low <- ci_up <- mat_temp <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
      cis <- apply(x$pcor_samples, 2, quantile, c(low, up))
      ci_low[] <- cis[1,]
      ci_up[] <- cis[2,]

      summ <- data.frame(edge = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         temp1 = ci_low[upper.tri(ci_low)],
                         temp2 = ci_up[upper.tri(ci_up)],
                         check.names = F)

      colnames(summ) <- c("Edge", "Estimate", "Est.Error",  paste(c("lb.", "ub."),
                           gsub("*0.","", formatC( round(x$ci, 4), format='f', digits=2)), "%", sep = ""))
      cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Estimates: \n \n")
      print(summ, row.names = F,...)
      cat("--- \n")

    }else{
      cat("Probability:", x$prob, "\n")
      cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Pr.out: post prob outside of rope \n")
      cat("Pr.in: post prob inside of rope \n")
      cat("--- \n")

      p <- ncol(x$partials_non_zero)
      mat_names <- mu_mat <- rope_in  <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]

      rope_in[] <- x$in_rope

      cat("Estimates: \n \n")
      summ <- data.frame(edge = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         "pr_out" = 1 - rope_in[upper.tri(rope_in)],
                         "pr_in" = rope_in[upper.tri(rope_in)],
                         check.names = F)

      colnames(summ) <- c("Edge", "Estimate",
                          "Est.Error",  "Pr.out", "Pr.in")
      print(summ, row.names = F,...)
      cat("--- \n")
    }
  }
}
