#' Title
#'
#' @param qData 
#' @param sTab 
#' @param VarRubin 
#' @param comp.type 
#' @param robust 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
mi4limma <- function (qData, sTab, VarRubin, comp.type = "OnevsOne", robust = F) 
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
