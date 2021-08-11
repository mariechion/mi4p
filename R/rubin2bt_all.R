# Computes the between-imputation component in the 2nd Rubin's rule
# for all peptides

#' Title
#'
#' @param data 
#' @param funcmean 
#' @param metacond 
#' @param is.parallel 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
rubin2bt.all <- function(data,funcmean = meanImp_emmeans, metacond, is.parallel = FALSE, verbose = FALSE){
  if (is.parallel){
    iforeach<-NULL
    requireNamespace("foreach",quietly = TRUE)
    res <- foreach::foreach(iforeach=1:dim(data)[1], .errorhandling = 'remove', .verbose = verbose) %dopar% 
      rubin2bt.one(iforeach,data=data,funcmean=funcmean,metacond=metacond)
  }
  else {
    res <- lapply(1:dim(data)[1],rubin2bt.one,data=data,funcmean=funcmean,metacond=metacond)
  }
  return(res)
}
