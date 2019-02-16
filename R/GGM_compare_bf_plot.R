#' Title
#'
#' @param X
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_bf_plot <- function(X, cutoff){



  melt_01 <- reshape::melt(X$BF_01)
  melt_10 <- reshape::melt(X$BF_10)

  if(class(X) != "GGM_compare_bf")stop("incorrect object class")

  melt_01 <- subset(melt_01, value > cutoff)
  melt_01$X1 <- as.factor(melt_01$X1)

  max_01 <- log(max(melt_01$value))



  melt_10 <- subset(melt_10, value > cutoff)
  melt_10$X1 <- as.factor(melt_10$X1)
  max_10 <- log(max(melt_10$value))




  plt_null <- ggplot(melt_01, aes(x = as.factor(X2),
                                  y = as.factor(X1),
                                  fill = log(value))) +
    geom_tile() +
    scale_fill_gradientn(colours = c("white", "yellow", "red"),
                         values = c(0, .1, 1),
                         limits = c(log(cutoff), max_01),
                         labels = round(c(3,   exp(max_01))),
                         breaks =  c(log(cutoff),  max_01),
                         name = "BF 01") +

    xlab("Null Hypothesis Matrix") + ylab(" ") +
    ylim(rev(levels(melt_01$X1)))+
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())  +
    theme(panel.grid = element_blank(),
          legend.title=element_text(size=9),
          legend.text=element_text(size=9))




  plt_alt <- ggplot(melt_10, aes(x = as.factor(X2),
                                 y = as.factor(X1),
                                 fill = round(log(value)))) +
    geom_tile() +
    scale_fill_gradientn(colours = c("white", "lightblue", "purple"),
                         values = c(0, .1, 1), limits = round(c(log(3), max_10)),
                         labels = round(c(3,   exp(max_10))),
                         breaks =  round(c(log(cutoff),  max_10)),
                         name = "BF 10") +


    xlab("Alternative Hypothesis Matrix") + ylab(" ") +
    ylim(rev(levels(melt_10$X1)))+
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.title=element_text(size=9),
          legend.text=element_text(size=9) )

  #plt_together <- cowplot::plot_grid(plt_null, plt_alt)


  #plt_together
  list(plt_null = plt_null, plt_alt = plt_alt)
}
