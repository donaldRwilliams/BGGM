#' Generate Ordinal and Binary data
#'
#' @param n number of observations
#' @param p number of variables
#' @param levels number of categories
#' @param cor_mat true correlation matrix
#' @param empirical logical. If true, mu and Sigma specify the empirical not
#'                  population mean and covariance matrix.
#'
#' @return
#' @export
#'
#' @examples
#'
#' main <- BGGM::ptsd_cor1[1:5,1:5]
#' p <- ncol(main)
#'
#' pcors <- -(cov2cor(solve(main)) -diag(p))
#' diag(pcors) <- 1
#' pcors <- ifelse(abs(pcors) < 0.05, 0, pcors)
#'
#' inv <-  -pcors
#' diag(inv) <- 1
#' cors <- cov2cor( solve(inv))
#'
#' # example data
#' Y <- BGGM::gen_ordinal(n = 500, p = 5,
#'                    levels = 2,
#'                     cor_mat = cors,
#'                     empirical = FALSE)

gen_ordinal <- function(n, p, levels,  cor_mat, empirical = FALSE){
  ls <- list()
  for(i in 1:p){
    temp <- table(sample(c(1:levels),
                         size = n,
                         replace = T))
    ls[[i]] <- as.numeric(temp / sum(temp))
  }
  junk <- capture.output(data <- rmvord_naiv(n = n, probs =  ls,
                                             Cor = cor_mat,
                                             empirical = empirical))
  data
}
