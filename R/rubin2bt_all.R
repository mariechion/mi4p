#' @title 2nd Rubin's rule Between-Imputation component (all peptides)
#' 
#' @description Computes the between-imputation component in the 2nd Rubin's rule for all peptides.
#'
#' @param data dataset
#' @param funcmean function that should be used to compute the mean
#' @param metacond a factor to specify the groups
#' @param is.parallel should parallel computing be used?
#' @param verbose should messages be displayed?
#'
#' @return List of variance-covariance matrices.
#' 
#' @author Frédéric Bertrand
#' 
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}. arxiv:2108.07086. \url{https://arxiv.org/abs/2108.07086}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin2bt.all(datasim_imp[1:5,,],funcmean = meanImp_emmeans,
#' attr(datasim,"metadata")$Condition)
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
