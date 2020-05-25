#' @title Partial Correlation Sum
#'
#' @name pcor_sum
#'
#' @description Compute and test partial correlation sums either within or between GGMs
#' (e.g., different groups), resulting in a posterior distribution.
#'
#' @param ... An object of class \code{estimate}. This can be either one or two fitted objects.
#'
#' @param relations Character string. Which partial correlations should be summed?
#'
#' @param iter Number of iterations (posterior samples; defaults to the number in the object).
#'
#' @return An object of class \code{posterior_sum}, including the sum and possibly the difference for
#' two sums.
#'
#' @details
#' Some care must be taken when writing the string for \code{partial_sum}. Below are several examples
#'
#' \strong{Just a Sum}:
#' Perhaps a sum is of interest, and not necessarily the difference of two sums. This can be written as
#'
#' \itemize{
#' \item \code{partial_sum <-  c("A1--A2 + A1--A3 + A1--A4")}
#' }
#'
#' which will sum those relations.
#'
#' \strong{Comparing Sums}:
#' When comparing sums, each must be seperated by "\code{;}". For example,
#'
#' \itemize{
#' \item \code{partial_sum <-  c("A1--A2 + A1--A3; A1--A2 + A1--A4")}
#' }
#'
#' which will sum both and compute the difference. Note that there cannot be more than two sums, such
#' that \code{c("A1--A2 + A1--A3; A1--A2 + A1--A4; A1--A2 + A1--A5")} will result in an error.
#'
#' \strong{Comparing Groups}:
#'
#' When more than one fitted object is suppled to \code{object} it is assumed that the groups
#' should be compared for the same sum. Hence, in this case, only the sum needs to be written.
#'
#' \itemize{
#' \item \code{partial_sum <-  c("A1--A2 + A1--A3 + A1--A4")}
#' }
#'
#' The above results in that sum being computed for each group and then compared.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- bfi
#'
#' # males
#' Y_males <- subset(Y, gender == 1, select = -c(education, gender))[,1:5]
#'
#' # females
#' Y_females <- subset(Y, gender == 2, select = -c(education, gender))[,1:5]
#'
#' # males
#' fit_males <- estimate(Y_males, seed = 1)
#'
#' # fit females
#' fit_females <- estimate(Y_females, seed = 2)
#'
#'
#' sums <- pcor_sum(fit_males,
#'                  fit_females,
#'                  relations = "A1--A2 + A1--A3")
#' # print
#' sums
#'
#' # plot difference
#' plot(sums)[[3]]
#' }
pcor_sum <- function(..., iter = NULL, relations){

  # collect ...
  collect_objects <- list(...)

  # number of groups
  groups <- length(collect_objects)

  # partial_sum_i
  partial_sum_i <- list()

  if(is.null(iter)){
    iter <- collect_objects[[1]]$iter
  }

  # separate to count
  count_sums <- strsplit(relations, "\\;")[[1]]

  # how many sums ?
  n_sums <- length(count_sums)

  remove_space <- gsub("[[:space:]]", "", count_sums)

  remove_plus <- gsub("\\+", replacement = " ",  remove_space)

  each_sum <- strsplit(remove_plus, split = "[[:space:]]")

  if (n_sums > 2) {
    stop("there is only support for 'at most' two sums")
  }

  # start one group
  if (groups == 1) {

    if(!all(c("estimate", "default") %in% class(collect_objects[[1]]))){
      stop("the object must be of class 'estimate'")
    }



    # posterior samples
    samps <- posterior_samples(collect_objects[[1]])[1:iter,]

    if (n_sums == 1) {

      # sum
      sums <- lapply(1:1, function(x) {

        sum_i <-  eval(parse(text = paste0("samps[,'",
                                           each_sum[[x]], "']",
                                           collapse = "+")))
      })

      # assign name for printing
      names(sums) <- remove_space
      diff <- NULL

      # start 2 sums
    } else {

      sums <- lapply(1:2, function(x) {
        sum_i <-  eval(parse(text = paste0(
          "samps[,'",
          each_sum[[x]], "']",
          collapse = "+"
        )))
      })

      diff <- sums[[1]] - sums[[2]]
      names(sums) <- remove_space

    }
  } else if (groups == 2) {

    if (!all(c("estimate", "default") %in% class(collect_objects[[1]]))) {
      stop("the object must be of class 'estimate'")
    }

    if (!all(c("estimate", "default") %in% class(collect_objects[[2]]))) {
      stop("the object must be of class 'estimate'")
    }

    if (n_sums > 1) {
      stop("only one sum can be specified when there are two groups")
    }

    sums <- lapply(1:2, function(g) {

      samps <- posterior_samples(collect_objects[[g]])[1:iter, ]

      sapply(1:1, function(x) {
        eval(parse(text = paste0(
          "samps[,'",
          each_sum[[x]], "']",
          collapse = "+"
        )))
      })
    })

    names(sums) <- paste0("g", 1:2, ": ", remove_space)
    diff <- sums[[1]] - sums[[2]]

  } else{
    stop("too many groups. only two is currently support")
  }

  partial_sum_i <- list(post_diff = diff,
                        post_sums = sums,
                        n_sums = n_sums,
                        iter = iter)



  returned_object <-  partial_sum_i
  class(returned_object) <- c("BGGM", "pcor_sum")
  return(returned_object)
}


print_pcor_sum <- function(x, cred = 0.95, row_names = TRUE){


  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Network Stats: Posterior Sum\n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("--- \n")
  cat("Estimates \n\n")

  # lower bound
  lb <- (1 - cred) / 2

  # upper bound
  ub <- 1 - lb

  if(is.null(x$post_diff)){
    cat("Sum:", "\n")
    res <- round(
      data.frame(Post.mean = mean(x$post_sums[[1]]),
                 Post.sd =  sd(x$post_sums[[1]]),
                 Cred.lb = quantile(x$post_sums[[1]], probs = lb),
                 Cred.ub = quantile(x$post_sums[[1]], probs = ub)
      ), 3)

    if(isTRUE(row_names)){
      rownames(res) <- names(x$post_sums)
    } else {
      rownames(res) <- NULL
    }
    print(res, row.names = row_names)
  } else {
    cat("Sum:", "\n")
    dat_i <- list()

    for(i in 1:2){
      dat_i[[i]] <-  round(
        data.frame(Post.mean = mean(x$post_sums[[i]]),
                   Post.sd =  sd(x$post_sums[[i]]),
                   Cred.lb = quantile(x$post_sums[[i]], probs = lb),
                   Cred.ub = quantile(x$post_sums[[i]], probs = ub)
        ), 3)

    }

    diff_res <- round(
      data.frame(Post.mean = mean(x$post_diff),
                 Post.sd =  sd(x$post_diff),
                 Cred.lb = quantile(x$post_diff, probs = lb),
                 Cred.ub = quantile(x$post_diff, probs = ub),
                 Prob.greater = mean(x$post_diff > 0),
                 Prob.less = mean(x$post_diff < 0)
      ), 3)
    res <- do.call(rbind.data.frame, dat_i)

    if(isTRUE(row_names)){

      rownames(res) <- names(x$post_sums)

    } else {
      rownames(res) <- NULL

    }
    rownames(diff_res) <- NULL
    print(res, row.names = row_names)
    cat("--- \n\n")
    cat("Difference:\n")
    cat(paste(names(x$post_sums)[1]), "-", paste(names(x$post_sums)[2]), "\n\n")
    print(diff_res, row.names = FALSE)
    cat("--- \n")
  }
}

#' @title Plot \code{pcor_sum} Object
#'
#' @name plot.pcor_sum
#'
#' @param x An object of class \code{posterior_sum}
#'
#' @param fill Character string. What fill for the histogram
#'        (defaults to colorblind "pink")?
#'
#' @param ... Currently ignored.
#'
#' @return A list of \code{ggplot} objects
#'
#' @export
#'
#' @note
#' \strong{Examples}:
#'
#' @seealso posterior_sum
plot.pcor_sum <- function(x,
                          fill = "#CC79A7",
                          ...){


  if(is.null( x$post_diff)){

    g1 <- ggplot(data.frame(x = x$post_sums[[1]]),
                 aes(x = x)) +
      geom_histogram(color = "white",
                     fill = fill) +
      xlab(names(x$post_sums)[1])

    if(length( x$post_sums) == 2){

      g2 <- ggplot(data.frame(x = x$post_sums[[2]]),
                   aes(x = x)) +
        geom_histogram(color = "white",
                       fill = fill) +
        xlab(names(x$post_sums)[2])

      list(g1 = g1, g2  = g2)

    } else {

      list(g1 = g1)

    }

  } else {

    g1 <- ggplot(data.frame(x = x$post_sums[[1]]),
                 aes(x = x)) +
      geom_histogram(color = "white",
                     fill = fill) +
      xlab(names(x$post_sums)[1])

    g2 <- ggplot(data.frame(x = x$post_sums[[2]]),
                 aes(x = x)) +
      geom_histogram(color = "white",
                     fill = fill) +
      xlab(names(x$post_sums)[2])

    diff <- ggplot(data.frame(x = x$post_diff),
                   aes(x = x)) +
      geom_histogram(color = "white",
                     fill = fill) +
      xlab("Difference")

    suppressWarnings( list(g1 = g1, g2 = g2, diff = diff))
  }


}

