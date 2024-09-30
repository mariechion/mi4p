#' @title 2nd Rubin's rule Within-Variance Component (all peptides)
#' 
#' @description Computes the within-variance component in the 2nd Rubin's rule for all peptides.
#'
#' @param data dataset
#' @param funcvar function that should be used to compute the variance
#' @param metacond a factor to specify the groups
#' @param is.parallel should parallel computing be used?
#' @param verbose should messages be displayed?
#'
#' @return List of variance-covariance matrices.
#' 
#' @author Frédéric Bertrand
#' 
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}.  \doi{doi:10.1371/journal.pcbi.1010420}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], 
#' conditions = attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin2wt.all(datasim_imp[1:5,,],funcvar = within_variance_comp_emmeans,
#' attr(datasim,"metadata")$Condition)
rubin2wt.all <- function(data, funcvar = mi4p::within_variance_comp_emmeans, 
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
  
