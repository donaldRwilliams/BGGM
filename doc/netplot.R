## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # need the developmental version
#  if (!requireNamespace("remotes")) {
#    install.packages("remotes")
#  }
#  
#  # install from github
#  remotes::install_github("donaldRwilliams/BGGM")
#  library(BGGM)
#  library(cowplot)

## ---- echo=FALSE, message=FALSE-----------------------------------------------
library(BGGM)

## ---- eval=FALSE--------------------------------------------------------------
#  # data
#  Y <- bfi[,1:25]
#  
#  # fit model
#  fit <- estimate(Y)

## ---- eval=FALSE--------------------------------------------------------------
#  # select the edge set
#  E <- select(fit,
#              cred = 0.95,
#              alternative = "two.sided")

## ---- eval=FALSE--------------------------------------------------------------
#  plot(E)

## ---- eval=FALSE--------------------------------------------------------------
#  plot(E,
#       # enlarge edges
#       edge_magnify = 5,
#       # cluster nodes
#       groups = comm,
#       # change layout
#       layout = "random")$plt +
#    # add custom labels
#    scale_color_brewer(breaks = c("A",
#                                  "C",
#                                  "E",
#                                  "N",
#                                  "O"),
#                       labels =   c("Agreeableness", "Conscientiousness",
#                                    "Extraversion", "Neuroticism",
#                                    "Opennness"),
#                       palette = "Set2")

