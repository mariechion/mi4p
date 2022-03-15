#' @title Computes a hierarchical differential analysis
#' 
#' @description Modified version of the \code{limmaCompleteTest} function from the \code{DAPAR} package to return both the fit and the results.
#' 
#' @param qData A matrix of quantitative data, without any missing values.
#' @param sTab A dataframe of experimental design (\code{pData()}).
#' @param comp.type A string that corresponds to the type of comparison. Values are: 'anova1way', 'OnevsOne' and 'OnevsAll'; default is 'OnevsOne'.
#'
#' @return A list of two dataframes : logFC and P_Value. The first one contains the logFC values of all the comparisons (one column for one comparison), the second one contains the pvalue of all the comparisons (one column for one comparison). The names of the columns for those two dataframes are identical and correspond to the description of the comparison.
#' @author Adapted from Hélène Borges, Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' @export
#'
#' @examples
#' set.seed(2016)
#' data(qData)
#' data(sTab)
#' limma <- limmaCompleteTest.mod(qData, sTab, comp.type='OnevsOne')
limmaCompleteTest.mod <- function (qData, sTab, comp.type = "OnevsOne") 
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
  design.matrix <- mi4p::make.design(sTab)
  if (!is.null(design.matrix)) {
    contra <- mi4p::make.contrast(design.matrix, condition = conds, 
                            contrast)
    cmtx <- limma::makeContrasts(contrasts = contra, levels = make.names(colnames(design.matrix)))
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(qData, 
                                                           design.matrix), cmtx))
    res.l <- mi4p::formatLimmaResult(fit, conds, contrast)
  }
  return(list(res.l = res.l, fit.s2 = fit$s2.post))
}
