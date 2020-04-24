params <-
list(EVAL = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
library(BGGM)

## -----------------------------------------------------------------------------
# data
Y <- BGGM::bfi[,1:5]

# sample posterior
fit <- estimate(Y, analytic = FALSE, 
                iter = 5000)

## -----------------------------------------------------------------------------
summary_fit <- summary(fit, cred = 0.90)

# print
summary_fit

## -----------------------------------------------------------------------------
summary_plot <- plot(summary_fit)

## -----------------------------------------------------------------------------
library(ggplot2)

# plot 
plot(summary_fit, width = 0) +
  # flip
  coord_flip() +
  # change the theme
  theme_bw() +
  # remove legend
  theme(legend.position = "none", 
        panel.grid  = element_blank()) +
  # add line at zero
  geom_hline(yintercept = 0, 
             linetype = "dotted", 
             alpha = 0.30) +
  # change color to black border
  geom_point(color = "black", size = 2) +
  # make inner point white
  geom_point(color = "white", size = 1.75) +
  # color blind
  scale_color_manual(values = c("#999999", "#CC79A7")) 

