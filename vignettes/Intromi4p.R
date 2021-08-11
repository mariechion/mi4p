## ----setup, include = FALSE---------------------------------------------------
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=FALSE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("mi4p")

## ---- eval = FALSE------------------------------------------------------------
#  devtools::install_github("mariechion/mi4p")

## -----------------------------------------------------------------------------
library(mi4p)

