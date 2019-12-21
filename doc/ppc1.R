params <-
list(EVAL = FALSE)

## ---- SETTINGS-knitr, include=FALSE-----------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  cache = TRUE,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.8,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
  )
library(BGGM)
library(ggplot2)

## ---------------------------------------------------------------------------------------
#  dat <- BGGM::bfi
#  
#  dat_males <- subset(dat, gender == 1)[,1:25]
#  
#  dat_female <- subset(dat, gender == 2)[,1:25]
#  
#  # fit model
#  fit1 <- ggm_compare_ppc(dat_males,
#                          dat_female,
#                          iter = 500,
#                          cores = 4)
#  

## ---------------------------------------------------------------------------------------
#  summary(fit1)

## ---------------------------------------------------------------------------------------
#  plot(fit1, critical = 0.05) +
#    theme_bw() +
#    theme(legend.position = "none")

