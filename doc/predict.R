params <-
list(EVAL = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
library(ggplot2)

## ---- message = FALSE---------------------------------------------------------
# packages
library(BGGM)
library(ggplot2)

# resilence data
# remove gender variable
dat <- subset(rsa, select = - gender)

## -----------------------------------------------------------------------------
# fit model
fit <- estimate(dat, iter = 1000)

## -----------------------------------------------------------------------------
# predict
pred <- predict(fit, summary = FALSE)

## -----------------------------------------------------------------------------
# summary
error <- mse(pred)

## -----------------------------------------------------------------------------
# plot 
plot(error)

## -----------------------------------------------------------------------------
plot(error) +
  theme_bw() +
  ggtitle("Predictability") +
  ylab("Mean Squared Error") +
  geom_point(size = 2, 
             color = "black") +
  geom_point(size = 1.5, 
             color = "white")

## ---- message=F---------------------------------------------------------------
fitted_pred <- plot(error, type = "ridgeline", 
     color = "red", 
     alpha =  0.75, 
     scale = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Node") +
  xlab("Mean Squared Error") +
  ggtitle("Predictability")
fitted_pred

## -----------------------------------------------------------------------------
pred <- posterior_predict(fit, iter = 250,
                          summary = FALSE)

## -----------------------------------------------------------------------------
error <- mse(pred)

# print summary
error

## -----------------------------------------------------------------------------
# plot 
plot(error)

## ---- message=F---------------------------------------------------------------
posterior_pred <- plot(error, type = "ridgeline", 
     color = "red", 
     alpha =  0.75, 
     scale = 2) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Node") +
  xlab("Mean Squared Error") +
  ggtitle("Predictability")
posterior_pred

## ---- warning=F, message=F----------------------------------------------------
top <- cowplot::plot_grid("", "", 
                          labels = c("Fitted", 
                              "Posterior Predictive"))

bottom <- cowplot::plot_grid(fitted_pred, 
                             posterior_pred)

cowplot::plot_grid(top, bottom, 
                   nrow = 2, 
                   rel_heights = c(1, 20))

