#' @title Differential analysis after multiple imputation
#'
#' @description 
#' This function performs hierarchical differential analysis using a moderated 
#' t-test statistic, which accounts for multiple imputation variability if 
#' applicable.
#'
#' @param qData A matrix of quantitative data, without any missing values. It 
#' should be the averaged matrix from the array resulting from 
#' \code{\link{multi.impute}}.
#' @param sTab The experimental matrix, also corresponding to the pData function 
#' of MSnbase.
#' @param VarRubin A numerical vector, resulting from \code{\link{proj_matrix}}.
#' It denotes the vector of projected variance-covariance matrices. It should be
#' of length the number of peptides or proteins considered.
#' @param comp.type A string that corresponds to the type of comparison. Values 
#' are: 'anova1way', 'OnevsOne' and 'OnevsAll'; default is 'OnevsOne'.
#' @param robust logical, should the estimation of df.prior and var.prior be 
#' robustified against outlier sample variances? (as in limma's eBayes)
#'
#'
#' @return A list of two dataframes : logFC and P_Value. The first one contains 
#' the logFC values of all the comparisons (one column for one comparison), the 
#' second one contains the pvalue of all the comparisons (one column for one 
#' comparison). The names of the columns for those two dataframes are identical 
#' and correspond to the description of the comparison.
#' 
#' @author 
#' Adapted by Marie Chion, from \code{limmaCompleteTest} of the 
#' \code{DAPAR} package by Hélène Borges, Thomas Burger, 
#' Quentin Giai-Gianetto and Samuel Wieczorek.
#' 
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting 
#' for multiple imputation-induced variability for differential analysis in 
#' mass spectrometry-based label-free quantitative proteomics}. 
#' arxiv:2108.07086. \url{https://arxiv.org/abs/2108.07086}.
#' 
#' @export
#'
#' @examples
#' set.seed(2016)
#' data(qData)
#' data(sTab)
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
  design.matrix <- mi4p::make.design(sTab)
  if (!is.null(design.matrix)) {
    contra <- mi4p::make.contrast(design.matrix, condition = conds, 
                            contrast)
    cmtx <- limma::makeContrasts(contrasts = contra, levels = make.names(colnames(design.matrix)))
    fit <- mi4p::eBayes.mod(limma::contrasts.fit(limma::lmFit(qData, 
                                                        design.matrix), cmtx), VarRubin, robust=robust)
    res.l <- mi4p::formatLimmaResult(fit, conds, contrast)
  }
  return(res.l)
}
