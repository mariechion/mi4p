## ----setup, include = FALSE---------------------------------------------------
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=FALSE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5
)

