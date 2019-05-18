#' Title
#'
#' @param x
#' @param hjust1
#' @param hjust2
#' @param hjust3
#' @param base_text
#' @param leg_text
#'
#' @return
#' @export
#'
#' @examples
plot.edge_compare.explore <- function(x, stack = TRUE, hjust1 = -1,
                                      hjust2 = -1, hjust3 = -1,
                                      base_text = 14,
                                      leg_text = 10, spread = 0.95){


  if(isFALSE(stack)){

    plts <- list()

    # for exhaustive
    if(x$alternative == "exhaustive"){

      for(i in 1:nrow(x$results)){

        # ith row
        dat_pie <- data.frame(hyp = colnames(x$results[i,4:6]),
                              prob = as.numeric(x$results[i,4:6]))
        # plot
        plts[[i]] <- ggplot(dat_pie, aes(y= prob,
                                         fill = hyp,
                                         x = hyp)) +
          # geom bar
          geom_bar(
            # size o bar
            width = spread,
            # plot probability
            stat = "identity",
            # unstack
            position = position_dodge(1)
          ) +
          # flip
          coord_flip() +
          # nice theme
          theme_classic(
            # text size
            base_size = base_text
          ) +
          # remove gaps
          scale_y_continuous(
            # remove gaps
            expand = c(0, 0),
            # span 0 to 1
            limits = c(0, 1)
          ) +
          # x label
          xlab("Hypothesis") +
          # legend as the hyp
          scale_fill_manual(
            name = "",
            # breaks in legend
            breaks = dat_pie$hyp,
            #
            labels = c(
              # Ho
              expression(
                H[0]:rho[diff]*" = "* 0),
              # H1
              expression(
                H[1]:rho[diff]*" > "* 0),
              # H2
              expression(
                H[2]:rho[diff]*" < "* 0)
            ),
            # colors
            values = c("#D55E00", "#CC79A7","#009E73" )) +
          # plot options
          theme(
            # top
            legend.position = "top",
            # right margin for the case of 100 %
            plot.margin = margin(0.1, .5, .1, .1, "cm"),
            # legend text size
            legend.text = element_text(size = leg_text)
          ) +
          # y label
          ylab(
            paste("Posterior Probability\n","(Contrast:",
                  x$results$contrast[i],")", sep = "")
          ) +
          # probs as text
          geom_text(aes(
            label = prob,
            y = prob,
            x = hyp),
            position = position_dodge(1),
            hjust = c(hjust1,hjust2, hjust3)
          )
      }
    }
    if(x$alternative == "greater"){



      for(i in 1:nrow(x$results)){

        # ith row
        dat_pie <- data.frame(hyp = colnames(x$results[i,5:6]),
                              prob = round(as.numeric(x$results[i,5:6]),3))
        # plot
        plts[[i]] <- ggplot(dat_pie, aes(y= prob,
                                         fill = hyp,
                                         x = hyp)) +
          # geom bar
          geom_bar(
            # size o bar
            width = spread,
            # plot probability
            stat = "identity",
            # unstack
            position = position_dodge(1)
          ) +
          # flip
          coord_flip() +
          # nice theme
          theme_classic(
            # text size
            base_size = base_text
          ) +
          # remove gaps
          scale_y_continuous(
            # remove gaps
            expand = c(0, 0),
            # span 0 to 1
            limits = c(0, 1)
          ) +
          # x label
          xlab("Hypothesis") +
          # legend as the hyp
          scale_fill_manual(
            name = "",
            # breaks in legend
            breaks = dat_pie$hyp,
            #
            labels = c(
              # Ho
              expression(
                H[0]:rho[diff]*" = "* 0),
              # H1
              expression(
                H[1]:rho[diff]*" > "* 0)
            ),
            # colors
            values = c("#D55E00", "#CC79A7" )) +
          # plot options
          theme(
            # top
            legend.position = "top",
            # right margin for the case of 100 %
            plot.margin = margin(0.1, .5, .1, .1, "cm"),
            # legend text size
            legend.text = element_text(size = leg_text)
          ) +
          # y label
          ylab(
            paste("Posterior Probability\n","(Contrast:",
                  x$results$contrast[i],")", sep = "")
          ) +
          # probs as text
          geom_text(aes(
            label = prob,
            y = prob,
            x = hyp),
            position = position_dodge(1),
            hjust = c(hjust1,hjust2)
          )
      }

    }

    if(x$alternative == "less"){

      for(i in 1:nrow(x$results)){

        # ith row
        dat_pie <- data.frame(hyp = colnames(x$results[i,5:6]),
                              prob = round(as.numeric(x$results[i,5:6]),3))
        # plot
        plts[[i]] <- ggplot(dat_pie, aes(y= prob,
                                         fill = hyp,
                                         x = hyp)) +
          # geom bar
          geom_bar(
            # size o bar
            width = spread,
            # plot probability
            stat = "identity",
            # unstack
            position = position_dodge(1)
          ) +
          # flip
          coord_flip() +
          # nice theme
          theme_classic(
            # text size
            base_size = base_text
          ) +
          # remove gaps
          scale_y_continuous(
            # remove gaps
            expand = c(0, 0),
            # span 0 to 1
            limits = c(0, 1)
          ) +
          # x label
          xlab("Hypothesis") +
          # legend as the hyp
          scale_fill_manual(
            name = "",
            # breaks in legend
            breaks = dat_pie$hyp,
            #
            labels = c(
              # Ho
              expression(
                H[0]:rho[diff]*" = "* 0),
              # H1
              expression(
                H[1]:rho[diff]*" < "* 0)
            ),
            # colors
            values = c("#D55E00", "#CC79A7" )) +
          # plot options
          theme(
            # top
            legend.position = "top",
            # right margin for the case of 100 %
            plot.margin = margin(0.1, .5, .1, .1, "cm"),
            # legend text size
            legend.text = element_text(size = leg_text)
          ) +
          # y label
          ylab(
            paste("Posterior Probability\n","(Contrast:",
                  x$results$contrast[i],")", sep = "")
          ) +
          # probs as text
          geom_text(aes(
            label = prob,
            y = prob,
            x = hyp),
            position = position_dodge(1),
            hjust = c(hjust1,hjust2)
          )
      }
    }

    if(x$alternative == "two.sided"){

      for(i in 1:nrow(x$results)){

        # ith row
        dat_pie <- data.frame(hyp = colnames(x$results[i,5:6]),
                              prob = round(as.numeric(x$results[i,5:6]),3))
        # plot
        plts[[i]] <- ggplot(dat_pie, aes(y= prob,
                                         fill = hyp,
                                         x = hyp)) +
          # geom bar
          geom_bar(
            # size o bar
            width = spread,
            # plot probability
            stat = "identity",
            # unstack
            position = position_dodge(1)
          ) +
          # flip
          coord_flip() +
          # nice theme
          theme_classic(
            # text size
            base_size = base_text
          ) +
          # remove gaps
          scale_y_continuous(
            # remove gaps
            expand = c(0, 0),
            # span 0 to 1
            limits = c(0, 1)
          ) +
          # x label
          xlab("Hypothesis") +
          # legend as the hyp
          scale_fill_manual(
            name = "",
            # breaks in legend
            breaks = dat_pie$hyp,
            #
            labels = c(
              # Ho
              expression(
                H[0]:rho[diff]*" = "* 0),
              # H1
              expression(
                H[1]:rho[diff]*" != "* 0)
            ),
            # colors
            values = c("#D55E00", "#CC79A7" )) +
          # plot options
          theme(
            # top
            legend.position = "top",
            # right margin for the case of 100 %
            plot.margin = margin(0.1, .5, .1, .1, "cm"),
            # legend text size
            legend.text = element_text(size = leg_text)
          ) +
          # y label
          ylab(
            paste("Posterior Probability\n","(Contrast:",
                  x$results$contrast[i],")", sep = "")
          ) +
          # probs as text
          geom_text(aes(
            label = prob,
            y = prob,
            x = hyp),
            position = position_dodge(1),
            hjust = c(hjust1,hjust2)
          )
      }

    }

    names(plts) <- x$results$contrast
  } else{

    if(x$alternative == "exhaustive") {

      probs <- unlist(lapply(1:length(x$results$contrast), function(z) unlist(x$results[z,4:6])))


      dat_pie <- data.frame(Contrast =rep(x$results$contrast, each = 3),
                            hyp = colnames(x$results[,4:6]),
                            prob = probs)


      # plot
      plts <- ggplot(dat_pie, aes(y= prob, fill = hyp,
                                  x = Contrast)) +

        ylab("Posterior Probability") +

        # geom bar
        geom_bar(
          # size o bar
          width = spread,
          # plot probability
          stat = "identity",
          # unstack
          # position = position_dodge(1)
        ) +
        coord_flip() +
        theme_classic(
          # text size
          base_size = base_text
        ) +
        # remove gaps
        scale_y_continuous(
          # remove gaps
          expand = c(0, 0),
          # span 0 to 1
          limits = c(0, 1)
        ) +
        theme(
          # top
          legend.position = "top",
          # right margin for the case of 100 %
          plot.margin = margin(0.1, .5, .1, .1, "cm"),
          # legend text size
          legend.text = element_text(size = leg_text)
        ) +
        scale_fill_manual(
          name = "",
          # breaks in legend
          breaks = dat_pie$hyp,
          #
          labels = c(
            # Ho
            expression(
              " "*H[0]:rho[diff]*" = "* 0),
            # H1
            expression(
              " "*H[1]:rho[diff]*" > "* 0),
            # H2
            expression(
              " "*H[2]:rho[diff]*" < "* 0)
          ),
          # colors
          values = c("#D55E00", "#CC79A7","#009E73" ))
    }

    if(x$alternative == "greater"){

      probs <- unlist(lapply(1:length(x$results$contrast), function(z) unlist(x$results[z,5:6])))

      dat_pie <- data.frame(Contrast =rep(x$results$contrast, each = 2),
                            hyp = colnames(x$results[,5:6]),
                            prob = probs)

      # plot
      plts <-  ggplot(dat_pie, aes(y= prob, fill = hyp,
                                   x = Contrast)) +

        ylab("Posterior Probability") +

        # geom bar
        geom_bar(
          # size o bar
          width = spread,
          # plot probability
          stat = "identity",
          # unstack
          # position = position_dodge(1)
        ) +
        coord_flip() +
        theme_classic(
          # text size
          base_size = base_text
        ) +
        # remove gaps
        scale_y_continuous(
          # remove gaps
          expand = c(0, 0),
          # span 0 to 1
          limits = c(0, 1)
        ) +
        theme(
          # top
          legend.position = "top",
          # right margin for the case of 100 %
          plot.margin = margin(0.1, .5, .1, .1, "cm"),
          # legend text size
          legend.text = element_text(size = leg_text)
        ) +
        scale_fill_manual(
          name = "",
          # breaks in legend
          breaks = dat_pie$hyp,
          #
          labels = c(
            # Ho
            expression(
              " "*H[0]:rho[diff]*" = "* 0),
            # H1
            expression(
              " "*H[1]:rho[diff]*" > "* 0)
          ),
          # colors
          values = c("#D55E00", "#CC79A7" ))



    }

    if(x$alternative == "less"){

      probs <- unlist(lapply(1:length(x$results$contrast), function(z) unlist(x$results[z,5:6])))


      dat_pie <- data.frame(Contrast =rep(x$results$contrast, each = 2),
                            hyp = colnames(x$results[,5:6]),
                            prob = probs)

      # plot
      plts <-  ggplot(dat_pie, aes(y= prob, fill = hyp,
                                   x = Contrast)) +

        ylab("Posterior Probability") +

        # geom bar
        geom_bar(
          # size o bar
          width = spread,
          # plot probability
          stat = "identity",
          # unstack
          # position = position_dodge(1)
        ) +
        coord_flip() +
        theme_classic(
          # text size
          base_size = base_text
        ) +
        # remove gaps
        scale_y_continuous(
          # remove gaps
          expand = c(0, 0),
          # span 0 to 1
          limits = c(0, 1)
        ) +
        theme(
          # top
          legend.position = "top",
          # right margin for the case of 100 %
          plot.margin = margin(0.1, .5, .1, .1, "cm"),
          # legend text size
          legend.text = element_text(size = leg_text)
        ) +
        scale_fill_manual(
          name = "",
          # breaks in legend
          breaks = dat_pie$hyp,
          #
          labels = c(
            # Ho
            expression(
              " "*H[0]:rho[diff]*" = "* 0),
            # H1
            expression(
              " "*H[1]:rho[diff]*" < "* 0)
          ),
          # colors
          values = c("#D55E00", "#CC79A7" ))
    }

    if(x$alternative == "two.sided"){

      probs <- unlist(lapply(1:length(x$results$contrast), function(z) unlist(x$results[z,5:6])))

      dat_pie <- data.frame(Contrast =rep(x$results$contrast, each = 2),
                            hyp = colnames(x$results[,5:6]),
                            prob = probs)


      # plot
      plts <-  ggplot(dat_pie, aes(y= prob, fill = hyp,
                                   x = Contrast)) +

        ylab("Posterior Probability") +

        # geom bar
        geom_bar(
          # size o bar
          width = spread,
          # plot probability
          stat = "identity",
          # unstack
          # position = position_dodge(1)
        ) +
        coord_flip() +
        theme_classic(
          # text size
          base_size = base_text
        ) +
        # remove gaps
        scale_y_continuous(
          # remove gaps
          expand = c(0, 0),
          # span 0 to 1
          limits = c(0, 1)
        ) +
        theme(
          # top
          legend.position = "top",
          # right margin for the case of 100 %
          plot.margin = margin(0.1, .5, .1, .1, "cm"),
          # legend text size
          legend.text = element_text(size = leg_text)
        ) +
        scale_fill_manual(
          name = "",
          # breaks in legend
          breaks = dat_pie$hyp,
          #
          labels = c(
            # Ho
            expression(
              " "*H[0]:rho[diff]*" = "* 0),
            # H1
            expression(
              " "*H[1]:rho[diff]*" != "* 0)
          ),
          # colors
          values = c("#D55E00", "#CC79A7" ))
    }
  }

  plts

}
