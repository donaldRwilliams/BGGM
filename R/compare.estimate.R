#' Edge (Partial COrrelation) Contrasts with Estimation Based Methods
#'
#' @description Allows for comparing partial correlations \emph{within} the same GGM--e.g., to determine the largest edge for each node. A difference can
#' either be asessed with a credible interval, or the region of practical equivalence. The latter allows for assessing practical equivalence--i.e.,
#' whether two edges have the same strength or partial correlation size.
#'
#' @param x object from \code{estimate} (\code{sampling = TRUE})
#' @param contrast partial correlations to compare
#' @param ci_width credible interval width
#' @param rope region of practical equivalence
#'
#' @return returned_object data.frame of results. use \code{summary} or \code{head}
#' @return call call used in \code{edge_difference}
#' @return ci credible interval used in \code{edge_difference}
#' @return rope region of practical equivalence used in \code{edge_difference}
#'
#' @return object of class \code{edge_compare.estimate}
#'
#' @note These contrasts are based on the posterior distribution, and credible intervals or the rope are used to determine differences. In the case of
#' the rope, it is also possible to assess similarity between edges. The function \code{edge_compare} is a generic, and for using the Bayes factor,
#' \code{x} should be an object from \code{explore}.
#'
#' @export
#'
#' @examples
#'
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample from posterior
#'
#' fit_sample <- estimate(Y, samples = 5000, analytic = F)
#'
#' edge_difference <- edge_compare(fit_sample,
#' contrast =  list("1--5 - 1--3",
#'                  "1--2 - 1--6",
#'                  "1--4 - 1--7",
#'                  "1--5 - 1--10",
#'                  "1--2 - 1--9"),
#'                  ci_width = 0.95,
#'                  rope = 0.1)
#'
#' head(edge_difference, nrow = 4)
edge_compare.estimate <- function(x, contrast, ci_width, rope = NULL){

  # if(!is.null(rope)){
  #   message("ci_width is ignored for decision rule, but used is for plotting")
  # }
  # lower interval
  low <- (1 - ci_width) / 2

  # upper interval
  up  <-  1 - low

  # number of variable
  p <- x$p

  # select partials
  pcors <- x$posterior_samples[, grep("pcors", colnames(x$posterior_samples))]

  # name partials to match contrasts
  colnames(pcors) <- unlist(lapply(1:p, function(x) paste(1:p, x, sep = "--")))

  options(warn=-1)
  if(length(contrast) == 1 && contrast != "all"){

    one <- sub(".* - ","", contrast)
    two <- gsub(" .*","", contrast)

    if(anyNA(match(c(one, two), colnames(pcors)))){
      stop("rewrite the contrasts. must be 1--2 - 1--3")
    }

    diff <- pcors[,two] - pcors[,one]


    # for rope (region of practical equivalence)
  if(is.numeric(rope)){

    # proportion in rope
    rp <- BGGM:::rope_helper(diff, rope)

    # returned object
    returned_object <-   data.frame( cbind(
                                # contrast name
                                "contrast" = contrast,
                                # posterior mean
                                post_mean = mean(diff),
                                # posterior sd
                                post_sd = sd(diff),
                                # credible interval
                                t(quantile(diff, c(low, up))),
                                # probability *out* rope
                                pr_out = 1 - rp,
                                # probability *in* rope
                                pr_in = rp), check.names = F)

    returned_object <- rapply(object = returned_object,
                              f = round,
                              classes = "numeric",
                              how = "replace",
                              digits = 4)
    }

  if(is.null(rope)){
  # no rope
  returned_object <- data.frame( cbind(
                              # contrast name
                              contrast = contrast,
                              # posterior mean
                              post_mean = mean(diff),
                              # posterior sd
                              post_sd = sd(diff),
                              # credible interval
                              t(quantile(diff, c(low, up)))), check.names  = F)


  returned_object <- rapply(object = returned_object,
                      f = round,
                      classes = "numeric",
                      how = "replace",
                      digits = 4)

  }

}
  # more than one contrast
  if(length(contrast) > 1){
    # store differences
    diff <- list()

    # store lists
    summ <- list()

    # check for list
    if(!is.list(contrast)){
      stop("more than one contrast must be a list")
    }
    for(i in 1:length(contrast)){

      # contrast


      one <- sub(".* - ","", contrast[i])
      two <- gsub(" .*","", contrast[i])


      if(anyNA(match(c(one, two), colnames(pcors)))){
        stop("rewrite the contrasts. must be 1--2 - 1--3")
      }


      # store differences
      diff[[i]] <- pcors[,two] - pcors[,one]

      # names
      names(diff)[[i]] <- contrast[i]

      if(is.numeric(rope)){
        # proportion in rope
        rp <- BGGM:::rope_helper(diff[[i]], rope)

        summ[[i]] <- data.frame( cbind(
                                      # contrast name
                                      "contrast" = contrast[i],
                                      # posterior mean
                                      post_mean = mean(diff[[i]]),
                                      # posterior sd
                                      post_sd = sd(diff[[i]]),
                                      # credible inteval
                                      t(quantile(diff[[i]], c(low, up))),
                                      # probability *out* rope
                                      pr_out = 1 - rp,
                                      # probability *in* rope
                                      pr_in = rp), check.names = F)

        summ[[i]] <- rapply(object = summ[[i]],
                            f = round,
                            classes = "numeric",
                            how = "replace",
                            digits = 4)

        } else{
      # no rope
      summ[[i]] <- data.frame(cbind(
                                    # contrast name
                                    "contrast" = contrast[i],
                                    # posterior mean
                                    post_mean = mean(diff[[i]]),
                                    # posterior sd
                                    post_sd = sd(diff[[i]]),
                                    # credible interval
                                    t(quantile(diff[[i]], c(low, up)))), check.names = F)

      summ[[i]] <- rapply(object = summ[[i]],
                          f = round,
                          classes = "numeric",
                          how = "replace",
                          digits = 4)
      }
    }
   # edge difference
   returned_object <-  do.call(rbind.data.frame, summ)


  }

  if(contrast == "all"){
    # store summaries
    summ <- list()

    # stored differences
    diff <- list()

    # pairwise combo's
    pairwise_comb <- t(combn(1:(.5 * (p *  (p - 1))), 2))

    # temporary matrix
    mat_names <- matrix(0, p , p)
    mat_names[] <- unlist( lapply(1:p, function(x) paste(1:p, x, sep = "--") ))

    # extract only upper triangular
    pcor_names <- mat_names[upper.tri(mat_names)]

    # for ith combinabtion
    for(i in 1:nrow(pairwise_comb)){

      # contrast
      con <- paste0(pcor_names[pairwise_comb[i,1]], "-", pcor_names[pairwise_comb[i,2]], collapse = "")

      # difference
      diff[[i]] <- pcors[,pcor_names[pairwise_comb[i,1]]] - pcors[,pcor_names[pairwise_comb[i,2]]]

      # names
      names(diff)[[i]] <- con

      if(is.numeric(rope)){
        # proportion in rope
        rp <- BGGM:::rope_helper(diff[[i]], rope = rope)


        summ[[i]] <- cbind.data.frame(
                                      # contrast
                                      contrast = con,
                                      # posterior mean
                                      post_mean = mean(diff[[i]]),
                                      # posterior sd
                                      post_sd = sd(diff[[i]]),
                                      # credible interval
                                      t(quantile(diff[[i]], c(low, up))),
                                      # probability *out* rope
                                      pr_out = 1- rp,
                                      # probability in rope
                                      pr_in = rp)

        summ[[i]] <- rapply(object = summ[[i]],
                            f = round,
                            classes = "numeric",
                            how = "replace",
                            digits = 3)

        }
      # no rope
      if(is.null(rope)){

        summ[[i]] <- cbind.data.frame(
                                      # contrast
                                      contrast = con,
                                      # posterior mean
                                      post_mean = mean(diff[[i]]),
                                      # posterior sd
                                      post_sd = sd(diff[[i]]),
                                      # credible interval
                                      t(quantile(diff[[i]], c(low, up))))

        summ[[i]] <- rapply(object = summ[[i]],
                            f = round,
                            classes = "numeric",
                            how = "replace",
                            digits = 3)

        }
    }

  returned_object <-    do.call(rbind.data.frame, summ)
}
  options(warn=0)
returned_object <- list(returned_object = returned_object,
                        call = match.call(),
                        ci = ci_width,
                        rope = rope, samples = diff)

class(returned_object) <- "edge_compare.estimate"

returned_object
}

