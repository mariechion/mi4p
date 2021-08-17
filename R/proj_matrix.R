#' @title Variance-Covariance Matrix Projection 
#' 
#' @description Use a projection of the given variance-covariance matrix.
#'
#' @param VarRubin.matrix A variance-covariance matrix.
#' @param metadata Metadata of the experiment.
#'
#' @return A list of variance-covariance matrices.
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}. arxiv:2108.07086. \url{https://arxiv.org/abs/2108.07086}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' VarRubin.matrix <- rubin2.all(datasim_imp[1:5,,],
#' attr(datasim,"metadata")$Condition)
#' proj_matrix(VarRubin.matrix, attr(datasim,"metadata"))
proj_matrix <- function(VarRubin.matrix, metadata) {
  return(as.numeric(lapply(VarRubin.matrix, function(aaa){
    DesMat = DAPAR::make.design(metadata)
    max(diag(aaa)%*%t(DesMat)%*%DesMat)
  })))
}