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
library(dplyr)
library(psych)
library(ggplot2)

## ---------------------------------------------------------------------------------------
#  library(psych)
#  dat <- bfi
#  
#  dat_males <- dat %>%
#    filter(gender == 1) %>%
#    select(1:25) %>%
#    na.omit()
#  
#  dat_female <- dat %>%
#    filter(gender == 2) %>%
#    select(1:25) %>%
#    na.omit()
#  
#  # fit model
#  fit1 <- ggm_compare_ppc(dat_males,
#                          dat_female,
#                          iter = 500,
#                          cores = 4)
#  

