# Computes first Rubin's rule for all peptides

#' Title
#'
#' @param data 
#' @param metacond 
#' @param funcmean 
#' @param is.parallel 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
rubin1.all <- function(data, metacond, funcmean = meanImp_emmeans, is.parallel = FALSE) {
  if (is.parallel) {
    iforeach<-NA
    requireNamespace("foreach",quietly = TRUE)
    res<-foreach::foreach(iforeach=1:dim(data)[1], .combine=cbind, 
                 .errorhandling = 'remove', .verbose = FALSE) %dopar% 
      rubin1.one(iforeach,data=data,
                 funcmean=funcmean,metacond=metacond)
    res<-t(simplify2array(res))
    rownames(res) <- rownames(data)
  }
  else {
    res <- t(simplify2array(lapply(1:dim(data)[1],
                                   rubin1.one,
                                   data=data,
                                   funcmean=funcmean,
                                   metacond=metacond)))
  }
  return(res)
}
