#' @title 2nd Rubin's rule Between-Imputation Component (a given peptide)
#'
#' @description Computes the between-imputation component in the 2nd Rubin's rule for a given peptide.
#' 
#' @param peptide peptide for which the variance-covariance matrix should be derived.
#' @param data dataset
#' @param funcmean function that should be used to compute the mean
#' @param metacond a factor to specify the groups
#'
#' @return A variance-covariance matrix.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' rubin2bt.one(1,datasim_imp,funcmean = meanImp_emmeans,
#' attr(datasim,"metadata")$Condition)
rubin2bt.one <- function(peptide,data,funcmean,metacond){
  funcmean_p = rubin1.one(peptide=peptide,data=data,funcmean = funcmean, metacond=metacond)
  outer_mat_prod = function(ind,peptide,funcmean_peptide){return(matrix((funcmean(ind,peptide=peptide,tabdata=data, metacond=metacond)-funcmean_peptide),ncol=1)%*%(funcmean(ind,peptide=peptide,tabdata=data, metacond=metacond)-funcmean_peptide))}
  return(1/(dim(data)[3]-1)*apply(simplify2array(lapply(1:dim(data)[3],outer_mat_prod,peptide=peptide,funcmean_p)),c(1,2),sum))    
}