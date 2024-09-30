#' @title 2nd Rubin's rule Within-Variance Component (a given peptide)
#'
#' @description Computes the within-variance component in the 2nd Rubin's rule for a given peptide.
#' 
#' @param peptide peptide for which the variance-covariance matrix should be derived.
#' @param data dataset
#' @param funcvar function that should be used to compute the variance
#' @param metacond a factor to specify the groups
#'
#' @return A variance-covariance matrix.
#' 
#' @author Frédéric Bertrand
#' 
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}.  \doi{doi:10.1371/journal.pcbi.1010420}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin2wt.one(1,datasim_imp,funcvar=within_variance_comp_emmeans,
#' attr(datasim,"metadata")$Condition)
rubin2wt.one <- function(peptide,data,funcvar,metacond){
  return(apply(simplify2array(lapply(1:dim(data)[3],
                                     funcvar,peptide=peptide,
                                     data=data,metacond=metacond)),c(1,2),mean))
}
