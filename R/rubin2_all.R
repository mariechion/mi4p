# Computes the 2nd Rubin's rule for all peptides

#' Title
#'
#' @param data 
#' @param metacond 
#' @param funcmean 
#' @param funcvar 
#' @param is.parallel 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
rubin2.all <- function(data, metacond, funcmean = meanImp_emmeans,
                       funcvar = within_variance_comp_emmeans, is.parallel = FALSE) {
  metacond = as.factor(metacond)
  if (is.null(dim(data)[3])) {
    data.array <- array(data = as.matrix(data), dim = c(nrow(data), ncol(data), 2))
  }
  else {
    data.array <- data
  }
  correct = (1+dim(data.array)[3])/dim(data.array)[3]
  Wp = rubin2wt.all(data = data.array, 
                             funcvar = funcvar, 
                             metacond = metacond, is.parallel)
  Bp = rubin2bt.all(data = data.array, 
                             funcmean = funcmean, 
                             metacond = metacond, is.parallel)
  if(is.parallel) {
    iforeach<-NULL
    requireNamespace("foreach",quietly = TRUE)
    res <- foreach::foreach(iforeach=1:length(Bp), .errorhandling = 'stop', .verbose = TRUE) %dopar% 
      (Wp[[iforeach]]+correct*Bp[[iforeach]])
  }
  else {
    res <- lapply(1:length(Bp), function(iforeach) {
      Wp[[iforeach]]+correct*Bp[[iforeach]]
    })
  }
  return(res)
}