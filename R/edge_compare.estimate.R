#' Title
#'
#' @param x
#' @param contrast
#' @param ci_width
#' @param rope
#'
#' @return
#' @export
#'
#' @examples
edge_compare.estimate <- function(x, contrast, ci_width, rope = NULL){

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

  if(length(contrast) == 1 & contrast != "all"){

    # remove spaces
    con <- gsub("[[:space:]]", "", contrast)

    # contrast 1
    one <- substr(con, 1, 4)

    # contrast 2
    two <- substr(con, 6, 9)

    # difference
    diff <- pcors[,one] - pcors[,two]


    # for rope (region of practical equivalence)
  if(is.numeric(rope)){

    # proportion in rope
    rp <- BGGM:::rope_helper(diff, rope)

    # returned object
    returned_object <- cbind.data.frame(
                                # contrast name
                                contrast = con,
                                # posterior mean
                                post_mean = mean(diff),
                                # posterior sd
                                post_sd = sd(diff),
                                # credible interval
                                t(quantile(diff, c(low, up))),
                                # probability *out* rope
                                pr_out_rope = 1 - rp,
                                # probability *in* rope
                                pr_n_rope = rp)

    returned_object <- rapply(object = returned_object,
                              f = round,
                              classes = "numeric",
                              how = "replace",
                              digits = 4)
    }

  if(is.null(rope)){
  # no rope
  returned_object <- cbind.data.frame(
                              # contrast name
                              contrast = con,
                              # posterior mean
                              post_mean = mean(diff),
                              # posterior sd
                              post_sd = sd(diff),
                              # credible interval
                              t(quantile(diff, c(low, up))))


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
      con <- gsub("[[:space:]]", "", contrast[i])

      # contrast 1
      one <- substr(con, 1, 4)

      # contrast 2
      two <- substr(con, 6, 9)

      # store differences
      diff[[i]] <- pcors[,one] - pcors[,two]

      # names
      names(diff)[[i]] <- con

      if(is.numeric(rope)){
        # proportion in rope
        rp <- BGGM:::rope_helper(diff[[i]], rope)

        summ[[i]] <- cbind.data.frame(
                                      # contrast name
                                      contrast = con,
                                      # posterior mean
                                      post_mean = mean(diff[[i]]),
                                      # posterior sd
                                      post_sd = sd(diff[[i]]),
                                      # credible inteval
                                      t(quantile(diff[[i]], c(low, up))),
                                      # probability *out* rope
                                      pr_out_rope = 1 - rp,
                                      # probability *in* rope
                                      pr_n_rope = rp)

        summ[[i]] <- rapply(object = summ[[i]],
                            f = round,
                            classes = "numeric",
                            how = "replace",
                            digits = 4)

        } else{
      # no rope
      summ[[i]] <- cbind.data.frame(
                                    # contrast name
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
                                      pr_out_rope = 1- rp,
                                      # probability in rope
                                      pr_n_rope = rp)

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

returned_object <- list(returned_object = returned_object,
                        call = match.call(),
                        ci = ci_width,
                        rope = rope, samples = diff)

class(returned_object) <- "edge_compare.estimate"

returned_object
}

