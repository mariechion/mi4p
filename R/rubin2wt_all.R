# Computes the within-variance component in the 2nd Rubin's rule 
# for all peptides

#' Title
#'
#' @param data 
#' @param funcvar 
#' @param metacond 
#' @param is.parallel 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
rubin2wt.all <- function(data, funcvar = within_variance_comp_emmeans, 
                         metacond, is.parallel = FALSE, verbose = TRUE) {
  if (is.parallel) {
    iforeach<-NULL
    requireNamespace("foreach",quietly = TRUE)
    res<-foreach::foreach(iforeach=1:dim(data)[1], .errorhandling = 'remove', .verbose = verbose) %dopar% 
      rubin2wt.one(iforeach,data=data,funcvar=funcvar,metacond=metacond)
  }
  else {
    res <- lapply(1:dim(data)[1],rubin2wt.one,data=data,funcvar=funcvar,metacond=metacond)
  }
  return(res)
}
  
