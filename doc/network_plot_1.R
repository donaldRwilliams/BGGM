params <-
list(EVAL = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
library(BGGM)

## -----------------------------------------------------------------------------
dat <- subset(tas, select = - gender)

fit <- explore(dat, 
               prior_sd = 0.25, 
               iter = 5000)

## -----------------------------------------------------------------------------
sel <- select(fit, BF_cut = 3)

## -----------------------------------------------------------------------------
plots <- plot(sel)

## -----------------------------------------------------------------------------
plots$plt

## -----------------------------------------------------------------------------
library(ggplot2)

# change node names 
# add labels (e.g., group the items)
# node names in white
# node text size = 6
# change node size
# increase edge width
# transparency
plots <- plot(sel, layout = "circle", 
              node_labels_color = "white",
              node_groups = BGGM:::tas_labels, 
              txt_size = 6, node_outer_size = 11, 
              node_inner_size = 8, 
              edge_multiplier = 5, 
              alpha = 0.3)
 # remove legend name and set palette
plots$plt + 
  scale_color_brewer(name = NULL, 
                     palette = "Set1")  +
    # add title
  ggtitle("Example Title") +
  # make title larger and add box around legend
  theme(plot.title = element_text(size = 20), 
        legend.background = element_rect(color = "black"))

## -----------------------------------------------------------------------------
# add different layout
plots <- plot(sel, "fruchtermanreingold", 
              node_labels_color = "white",
              node_groups = BGGM:::tas_labels, 
              txt_size = 6, node_outer_size = 11, 
              node_inner_size = 8, 
              edge_multiplier = 5, 
              alpha = 0.3)

# further customize
plots$plt + 
  # remove legend name and set palette
  scale_color_brewer(name = NULL, 
                     palette = "Set1")  +
  # add title
  ggtitle("Example Title") +
  # make title larger and add box around legend
  theme(plot.title = element_text(size = 20), 
        legend.background = element_rect(color = "black"))

## -----------------------------------------------------------------------------
plots$plt_null

