#' Plotting (node-wise) Differences in Predictive Accuracy
#'

#' @param x object of type \code{predict}
#' @param limits upper and lower boundary for the gradient
#'
#' @return
#' @export
#'
#' @examples
plot.compare.predict  <- function(x, limits){

  # check limits
  if(length(limits) != 2){
    stop("limites must be of length two--e.g., c(-1, 1)")
  }

  # extract 1st node in the pairwise comparsion
  one <- unlist(noquote(lapply(1:length(x$summary_error$contrast),
                function(y) sub(".*- ", "", x$summary_error$contrast[y]))))

  # extract 2nd node in the pairwise comparsion
  two <- unlist(noquote(lapply(1:length(x$summary_error$contrast),
                                function(y) gsub(" .*$", "", x$summary_error$contrast[y]))))

  # data.frame for plotting
  dat_plot <- rbind.data.frame(cbind.data.frame(
                                               # 1st node
                                               X1 = one,
                                               # 2nd node
                                               X2 = two,
                                               # R2 difference
                                               value = x$summary_error$post_mean) ,

                              cbind.data.frame(
                                               # 1st node
                                               X1 = two,
                                               # 2nd node
                                               X2 = one,
                                               # R2 difference
                                               value = x$summary_error$post_mean))

  # check interval for 0
  dat_plot$sig <- ifelse(x$summary_error[,4] < 0 & x$summary_error[,5] > 0, 0, 1)

  # create factor for ordering axes
  dat_plot$X1 <- factor(
                        dat_plot$X1,
                        levels = 1:length(levels(dat_plot$X1)),
                        labels = 1:length(levels(dat_plot$X1))
                        )

  # create factor for ordering axes
  dat_plot$X2 <- factor(dat_plot$X2,
                        levels = 1:length(levels(dat_plot$X1)),
                        labels = 1:length(levels(dat_plot$X1)))

  # predictive difference plot
  plt1 <- ggplot(dat_plot, aes(x = as.factor(X2),
                               y = as.factor(X1),
                               fill = value)) +
    # heat map
    geom_tile() +
    # gradient for R2 difference
    scale_fill_gradientn(colours = c("brown3", "white",  "palegreen3"),
                         name = "Difference",
                         # adjust legend
                         limits = limits) +
    # text size
    theme_bw(base_size = 12) +
    # remove grid
    theme(panel.grid = element_blank()) +
    # reverse y levels
    scale_y_discrete(limits= rev(levels(dat_plot$X1)),
                     expand = c(0, 0)) +
    # remove gaps around plot
    scale_x_discrete(expand = c(0,0)) +
    ylab("") +
    xlab("")

 # exclusion of 0 plot
 plt2 <- ggplot(dat_plot, aes(x = as.factor(X2),
                               y = as.factor(X1),
                               fill = as.factor(sig))) +
    # heat map
    geom_tile(show.legend = F) +
    # fill colors
    scale_fill_manual(values = c("white", "grey35")) +
    # text size
    theme_bw(base_size = 12) +
    # remove grid
    theme(panel.grid = element_blank()) +
    # reverse y levels
    scale_y_discrete(limits= rev(levels(dat_plot$X1)),
                     expand = c(0, 0)) +
    # remove gaps around plot
    scale_x_discrete(expand = c(0,0))     +
    ylab("") +
    xlab("")

  # plots
  plots <-  list(pred_diff = plt1, sig_diff = plt2)
  return(plots)
}






