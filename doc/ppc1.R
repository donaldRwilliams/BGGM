params <-
list(EVAL = FALSE)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 7, fig.align = "center")
library(ggplot2)
library(BGGM)

## -----------------------------------------------------------------------------
dat <- BGGM::bfi 

dat_males <- subset(dat, gender == 1)[,1:25]

dat_female <- subset(dat, gender == 2)[,1:25]

# fit model
fit1 <- ggm_compare_ppc(dat_males, 
                        dat_female, 
                        iter = 500, 
                        cores = 4)


## ----fit1---------------------------------------------------------------------
summary(fit1)

## -----------------------------------------------------------------------------
plot(fit1, critical = 0.05) + 
  theme_bw() + 
  theme(legend.position = "none")

