#' Title
#'
#' @param qData A matrix of quantitative data, without any missing values.
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' @param VarRubin A variance-covariance matrix
#' @param comp.type A string that corresponds to the type of comparison. Values are: 'anova1way', 'OnevsOne' and 'OnevsAll'; default is 'OnevsOne'.
#' @param robust logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances? (as in limma's eBayes)
#'
#' @return A list of two dataframes : logFC and P_Value. The first one contains the logFC values of all the comparisons (one column for one comparison), the second one contains the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes are identical and correspond to the description of the comparison.
#' @export
#'
#' @examples
#' #' library(DAPAR)
#' set.seed(2016)
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' qData <- Biobase::exprs(obj)
#' sTab <- Biobase::pData(obj)
#' fit.limma <- mi4limma(qData, sTab, diag(1,2))
mi4limma <- function (qData, sTab, VarRubin, comp.type = "OnevsOne", robust = FALSE) 
{
  switch(comp.type, OnevsOne = contrast <- 1, OnevsAll = contrast <- 2)
  rownames(sTab) <- sTab$Sample.name
  colnames(qData) <-sTab$Sample.name
  sTab.old <- sTab
  conds <- factor(sTab$Condition, levels = unique(sTab$Condition))
  sTab <- sTab[unlist(lapply(split(sTab, conds), function(x) {
    x["Sample.name"]
  })), ]
  qData <- qData[, unlist(lapply(split(sTab.old, conds), function(x) {
    x["Sample.name"]
  }))]
  conds <- conds[order(conds)]
  res.l <- NULL
  design.matrix <- DAPAR::make.design(sTab)
  if (!is.null(design.matrix)) {
    contra <- DAPAR::make.contrast(design.matrix, condition = conds, 
                            contrast)
    cmtx <- limma::makeContrasts(contrasts = contra, levels = make.names(colnames(design.matrix)))
    fit <- eBayes.mod(limma::contrasts.fit(limma::lmFit(qData, 
                                                        design.matrix), cmtx), VarRubin, robust=robust)
    res.l <- mi4p::formatLimmaResult(fit, conds, contrast)
  }
  return(res.l)
}
