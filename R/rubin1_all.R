#' @title First Rubin rule (all peptide)
#' 
#' @description Computes the first Rubin's rule for all the peptide.
#'
#' @param data dataset
#' @param metacond a factor to specify the groups
#' @param funcmean function that should be used to compute the mean
#' @param is.parallel should parallel computing be used?
#' @param verbose should messages be displayed?
#'
#' @return A vector of estimated parameters.
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}. arxiv:2108.07086. \url{https://arxiv.org/abs/2108.07086}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin1.all(datasim_imp[1:5,,],funcmean = meanImp_emmeans,
#' attr(datasim,"metadata")$Condition)
rubin1.all <- function(data, metacond, funcmean = meanImp_emmeans, is.parallel = FALSE, verbose=FALSE) {
  if (is.parallel) {
    iforeach<-NA
    requireNamespace("foreach",quietly = TRUE)
    res<-foreach::foreach(iforeach=1:dim(data)[1], .combine=cbind, 
                 .errorhandling = 'remove', .verbose = verbose) %dopar% 
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
