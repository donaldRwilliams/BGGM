#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.select.estimate <- function(x){
  if(is.null(x$rope)){
    partials <- reshape::melt(x$partials)
    partials$X1 <- as.factor(partials$X1)
    partials <- subset(partials, value != 0)


    plt <- ggplot(partials, aes(x = as.factor(X2),
                                y = as.factor(X1),
                                fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"), name =
                             "") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank()) +
      scale_y_discrete(limits= rev(levels(partials$X1)), expand = c(0, 0)) +
      scale_x_discrete(expand = c(0,0))     +
      ylab("") +
      xlab("")
  }


  if(!is.null(x$rope)){

    partials <- reshape::melt(x$partials)
    partials$X1 <- as.factor(partials$X1)
    partials <- subset(partials, value != 0)


    plt1 <- ggplot(partials, aes(x = as.factor(X2),
                                 y = as.factor(X1),
                                 fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"), name =
                             "") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(), plot.background = element_rect(fill  = "white")) +
      scale_y_discrete(limits= rev(levels(partials$X1)), expand = c(0, 0)) +
      scale_x_discrete(expand = c(0,0))     +
      ylab("") +
      xlab("")



    zeros <- reshape::melt(x$adjacency_zero)
    zeros$X1 <- as.factor(zeros$X1)
    # zeros <- subset(zeros, value != 0)



    plt2 <- ggplot(zeros, aes(x = as.factor(X2),
                              y = as.factor(X1),
                              fill = as.factor(value))) +
      geom_tile(show.legend = F) +
      scale_fill_manual(values = c("white", "grey35")) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank()) +
      scale_y_discrete(limits= rev(levels(partials$X1)), expand = c(0, 0)) +
      scale_x_discrete(expand = c(0,0))     +
      ylab("") +
      xlab("")


    plt <- list(plot_nonzero = plt1, plot_zero = plt2)







  }

  plt





}
