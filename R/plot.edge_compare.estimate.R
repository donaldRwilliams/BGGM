#' Title
#'
#' @param x
#' @param size
#' @param color
#' @param width
#' @param rope
#' @param prob
#'
#' @return
#' @export
#'
#' @examples
plot.edge_compare.estimate <- function(x,  size = 2, color = "red", width = .01, rope = NULL, prob = NULL){
  x <- x[1:4]
  if (is.null(rope)) {
    dat_plt <- x$returned_object
    dat_plt$sig <- ifelse(dat_plt$`2.5%` < 0 & dat_plt$`97.5%` > 0, 0, 1)
    dat_plt <- subset(dat_plt, sig == 1)
    dat_plt$contrast <- factor(dat_plt$contrast,
                               levels = dat_plt$contrast[order(dat_plt$post_mean)],
                               labels = dat_plt$contrast[order(dat_plt$post_mean)])

    plt <- ggplot(dat_plt, aes(y = post_mean,
                               x = contrast)) +
      geom_hline(yintercept = 0,
                 size = 1,
                 linetype = "dotted",
                 color = "grey50") +
      geom_point() +
      geom_errorbar(aes(ymin =  dat_plt$`2.5%`,
                        ymax = dat_plt$`97.5%`),
                    width = width) +
      geom_point(size = size,
                 color = color) +
      coord_flip() +
      xlab("Node") +
      ylab("Edge Difference") +
      theme_classic()
  }

  if(isTRUE(rope)){
    dat_plt <- x$returned_object

    if( sum( dat_plt$pr_out_rope > prob) == 0 &   sum( dat_plt$pr_n_rope  > prob) == 0 ){
      stop("no intervals exclude rope")
    }


    if(sum( dat_plt$pr_out_rope > prob) != 0 & sum( dat_plt$pr_n_rope > prob) != 0){

      dat_nonzero <- subset(dat_plt, dat_plt$pr_out_rope >= prob)
      dat_nonzero$contrast <- factor(dat_nonzero$contrast,
                                     levels = dat_nonzero$contrast[order(dat_nonzero$post_mean)],
                                     labels = dat_nonzero$contrast[order(dat_nonzero$post_mean)])


      dat_zero <- subset(dat_plt, dat_plt$pr_n_rope >= prob)
      dat_zero$contrast <- factor(dat_zero$contrast,
                                  levels = dat_zero$contrast[order(dat_zero$post_mean)],
                                  labels = dat_zero$contrast[order(dat_zero$post_mean)])


      plt_nonzero <- ggplot(dat_nonzero, aes(y = post_mean,

                                             x = contrast)) +

        geom_hline(yintercept = 0,
                   size = 1,
                   linetype = "dotted",
                   color = "grey50") +
        geom_point() +
        geom_errorbar(aes(ymin =  dat_nonzero$`2.5%`,
                          ymax = dat_nonzero$`97.5%`),
                      width = width) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Edge Difference") +
        theme_classic()  +
        annotate(geom = "rect", xmin = -Inf , xmax = Inf, ymax = -x$rope, ymin = x$rope, alpha = .1)

      plt_zero <- ggplot(dat_zero, aes(y = post_mean,

                                       x = contrast)) +

        geom_hline(yintercept = 0,
                   size = 1,
                   linetype = "dotted",
                   color = "grey50") +
        geom_point() +
        geom_errorbar(aes(ymin =  dat_zero$`2.5%`,
                          ymax = dat_zero$`97.5%`),
                      width = width) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Edge Difference") +
        theme_classic()  +
        annotate(geom = "rect", xmin = -Inf , xmax = Inf, ymax = -x$rope, ymin = x$rope, alpha = .1)


      plt <- list(plt_zero = plt_zero, plt_nonzero = plt_nonzero)

    }


    if(sum( dat_plt$pr_out_rope > prob) != 0){

      dat_nonzero <- subset(dat_plt, dat_plt$pr_out_rope >= prob)
      dat_nonzero$contrast <- factor(dat_nonzero$contrast,
                                     levels = dat_nonzero$contrast[order(dat_nonzero$post_mean)],
                                     labels = dat_nonzero$contrast[order(dat_nonzero$post_mean)])



      plt<- ggplot(dat_nonzero, aes(y = post_mean,

                                    x = contrast)) +

        geom_hline(yintercept = 0,
                   size = 1,
                   linetype = "dotted",
                   color = "grey50") +
        geom_point() +
        geom_errorbar(aes(ymin =  dat_nonzero$`2.5%`,
                          ymax = dat_nonzero$`97.5%`),
                      width = width) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Edge Difference") +
        theme_classic()  +
        annotate(geom = "rect",
                 xmin = -Inf,
                 xmax = Inf,
                 ymax = -x$rope,
                 ymin = x$rope,
                 alpha = .1)


    }

    if(sum( dat_plt$pr_n_rope > prob) != 0){

      plt <- ggplot(dat_zero, aes(y = post_mean,

                                  x = contrast)) +

        geom_hline(yintercept = 0,
                   size = 1,
                   linetype = "dotted",
                   color = "grey50") +
        geom_point() +
        geom_errorbar(aes(ymin =  dat_zero$`2.5%`,
                          ymax = dat_zero$`97.5%`),
                      width = width) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Edge Difference") +
        theme_classic()  +
        annotate(geom = "rect", xmin = -Inf , xmax = Inf, ymax = -x$rope, ymin = x$rope, alpha = .1)



    }

  }
  plt

}
