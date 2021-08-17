#' @title Multiple Imputation Estimate
#' 
#' @description Computes the multiple imputation parameter estimate using the emmeans package.
#'
#' @param ind index 
#' @param peptide name of the peptide
#' @param tabdata dataset
#' @param metacond a factor to specify the groups
#'
#' @return A vector.
#' @references M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple imputation-induced variability for differential analysis in mass spectrometry-based label-free quantitative proteomics}. arxiv:2108.07086. \url{https://arxiv.org/abs/2108.07086}.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' meanImp_emmeans(1,1,datasim_imp,attr(datasim,"metadata")$Condition)
meanImp_emmeans <- function(ind,peptide=1,tabdata,metacond) {
  tempdata <- cbind(reponse=tabdata[peptide,,ind],data.frame(fact=factor(metacond)))
  templm <- lm(reponse~fact,data=tempdata)
  return(summary(emmeans::emmeans(templm,specs = "fact"))$emmean)
}