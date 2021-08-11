# Modified version of limmaCompleteTest function from DAPAR package
# to return fit and res.

#' Title
#'
#' @param qData 
#' @param sTab 
#' @param comp.type 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
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
  design.matrix <- DAPAR::make.design(sTab)
  if (!is.null(design.matrix)) {
    contra <- DAPAR::make.contrast(design.matrix, condition = conds, 
                            contrast)
    cmtx <- limma::makeContrasts(contrasts = contra, levels = make.names(colnames(design.matrix)))
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(qData, 
                                                           design.matrix), cmtx))
    res.l <- mi4p::formatLimmaResult(fit, conds, contrast)
  }
  return(list(res.l = res.l, fit.s2 = fit$s2.post))
}
