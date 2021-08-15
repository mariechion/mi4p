#' @title Multiple Imputation Within Variance Component
#' 
#' @description Computes the multiple imputation within variance component using the emmeans package.
#'
#' @param ind index 
#' @param peptide name of the peptide
#' @param data dataset
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
#' within_variance_comp_emmeans(1,1,datasim_imp,
#' attr(datasim,"metadata")$Condition)
within_variance_comp_emmeans = function(ind,peptide,data,metacond){
  tempdata <- cbind(reponse=data[peptide,,ind],data.frame(fact=factor(metacond)))
  templm <- lm(reponse~fact,data=tempdata)
  return(vcov((emmeans::emmeans(templm,specs = "fact"))))
}
