params <-
list(EVAL = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
library(ggplot2)

## -----------------------------------------------------------------------------
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
pred <- predict(fit,  cred = 0.95)

## -----------------------------------------------------------------------------
# summary
print(pred)

## -----------------------------------------------------------------------------
# plot 
plot(pred)

## ---- message=FALSE, warning=FALSE--------------------------------------------
plot(pred) + 
  xlab("Item") +
  geom_point(aes(color = BGGM:::rsa_labels), 
             size = 4) +
   geom_point(size = 3, 
              color = "white") +
  scale_color_brewer(name = "Domain", 
                     palette = "Set1") +
  ggtitle("Predictability") +
  coord_cartesian() +
  theme(axis.text.x = element_text(angle = 60, 
                                   hjust = .5, 
                                   vjust = .5))

## -----------------------------------------------------------------------------
# training data
train_dat <- dat[1:600,]

# testing data
test_dat <- dat[601:675,]

# fit model
fit <- estimate(train_dat, iter = 1000)

## -----------------------------------------------------------------------------
pred <- predict(fit, 
                test_data = test_dat, 
                cred = 0.95)

pred

## -----------------------------------------------------------------------------
plot(pred) + 
  ggtitle("Out-of-Sample Predictability")

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(ggridges)
dat_res <- reshape2::melt(pred$post_samples)

dat_res$L1 <- factor(dat_res$L1, 
                     labels = order(pred$summary_error$post_mean, decreasing = F), 
                     levels = colnames(dat)[order(pred$summary_error$post_mean, decreasing = F)])

ggplot(dat_res, aes(x = value, y = L1)) + 
    geom_density_ridges(rel_min_height = 0.01) +
    theme_bw() +
    ylab("Item") +
    xlab("Bayesian R2") +
    theme(panel.grid.minor = element_blank())

