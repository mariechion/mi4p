#' @title Computes the 2nd Rubin's rule (all peptides)
#' 
#' @description Computes the total variance-covariance component in the 2nd Rubin's rule for all peptides.
#'
#' @param data dataset
#' @param metacond a factor to specify the groups
#' @param funcmean function that should be used to compute the mean
#' @param funcvar function that should be used to compute the variance
#' @param is.parallel should parallel computing be used?
#' @param verbose should messages be displayed?
#'
#' @return List of variance-covariance matrices.
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}. arxiv:2108.07086. \url{https://arxiv.org/abs/2108.07086}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin2.all(datasim_imp[1:5,,],attr(datasim,"metadata")$Condition)
rubin2.all <- function(data, metacond, funcmean = meanImp_emmeans,
                       funcvar = within_variance_comp_emmeans, is.parallel = FALSE, verbose = FALSE) {
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
                             metacond = metacond, is.parallel, verbose = verbose)
  Bp = rubin2bt.all(data = data.array, 
                             funcmean = funcmean, 
                             metacond = metacond, is.parallel, verbose = verbose)
  if(is.parallel) {
    iforeach<-NULL
    requireNamespace("foreach",quietly = TRUE)
    res <- foreach::foreach(iforeach=1:length(Bp), .errorhandling = 'stop', .verbose = verbose) %dopar% 
      (Wp[[iforeach]]+correct*Bp[[iforeach]])
  }
  else {
    res <- lapply(1:length(Bp), function(iforeach) {
      Wp[[iforeach]]+correct*Bp[[iforeach]]
    })
  }
  return(res)
}