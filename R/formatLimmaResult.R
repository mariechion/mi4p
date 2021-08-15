#' @title Format a Result from Limma
#' 
#' @description It is not exported by DAPAR and has to be reproduced here.
#' 
#' @param fit Limma fit
#' @param conds Condition vector
#' @param contrast Contrast vector
#'
#' @return A list of two dataframes : logFC and P_Value. The first one contains the logFC values of all the comparisons (one column for one comparison), the second one contains the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes are identical and correspond to the description of the comparison.
#' @author Adapted from the code of Samuel Wieczorek in the DAPAR package as it is an object that is not exported by the DAPAR package.
#' @export
#'
#' @examples
#' library(DAPAR)
#' set.seed(2016)
#' contrast=1
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' qData <- Biobase::exprs(obj)[1:1000,]
#' sTab <- Biobase::pData(obj)
#' rownames(sTab) <- sTab$Sample.name
#' colnames(qData) <-sTab$Sample.name
#' sTab.old <- sTab
#' conds <- factor(sTab$Condition, levels = unique(sTab$Condition))
#' sTab <- sTab[unlist(lapply(split(sTab, conds), function(x) {
#'   x["Sample.name"]
#' })), ]
#' qData <- qData[, unlist(lapply(split(sTab.old, conds), function(x) {
#'   x["Sample.name"]
#' }))]
#' conds <- conds[order(conds)]
#' res.l <- NULL
#' design.matrix <- DAPAR::make.design(sTab)
#' contra <- DAPAR::make.contrast(design.matrix, condition = conds, 
#'                                  contrast)
#' cmtx <- limma::makeContrasts(contrasts = contra, levels = make.names(colnames(design.matrix)))
#' fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(qData, 
#'                                                       design.matrix), cmtx))
#' res.l <- mi4p::formatLimmaResult(fit, conds, contrast)
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
