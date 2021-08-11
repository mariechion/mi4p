
#' Title
#' 
#' @description It is not exported by DAPAR and has to be reproduced here.
#' 
#' @param fit 
#' @param conds 
#' @param contrast 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
formatLimmaResult <- function (fit, conds, contrast) 
{
  res <- cbind(fit$coefficients, fit$p.value)
  Compa.Nb <- dim(fit$p.value)[2]
  cn <- c()
  for (i in 1:Compa.Nb) {
    if (contrast == 1) {
      compa <- stringr::str_match_all(colnames(fit$p.value)[i], 
                                      "[[:space:]]Condition([[:digit:]]+)")[[1]]
      cn[i] <- paste(unique(conds)[as.numeric(compa[1, 
                                                    2])], "_vs_", unique(conds)[as.numeric(compa[2, 
                                                                                                 2])], sep = "")
    }
    if (contrast == 2) {
      compa <- stringr::str_match_all(colnames(fit$p.value)[i], 
                                      "[[:space:]]Condition([[:digit:]]+)")[[1]]
      cn[i] <- paste(unique(conds)[as.numeric(compa[1, 
                                                    2])], "_vs_(all-", unique(conds)[as.numeric(compa[1, 
                                                                                                      2])], ")", sep = "")
    }
  }
  res.l <- list(logFC = as.data.frame(res[, 1:Compa.Nb]), P_Value = as.data.frame(res[, 
                                                                                      -(1:Compa.Nb)]))
  colnames(res.l$logFC) <- paste(cn, "logFC", sep = "_")
  colnames(res.l$P_Value) <- paste(cn, "pval", sep = "_")
  return(res.l)
}