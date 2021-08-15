#' @title First Rubin rule (a given peptide)
#' 
#' @description Computes the first Rubin's rule for a given peptide.
#'
#' @param peptide peptide for which the variance-covariance matrix should be derived.
#' @param data dataset
#' @param funcmean function that should be used to compute the mean
#' @param metacond a factor to specify the groups
#'
#' @return A vector of estimated parameters.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin1.one(1,datasim_imp,funcmean = meanImp_emmeans,
#' attr(datasim,"metadata")$Condition)
rubin1.one <- function(peptide,data,funcmean = meanImp_emmeans,metacond) {
  return(rowMeans(simplify2array(lapply(1:dim(data)[3],
                                        funcmean,
                                        tabdata=data,
                                        peptide=peptide,
                                        metacond=metacond))))
}