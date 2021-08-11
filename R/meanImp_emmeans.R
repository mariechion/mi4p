#' Title
#'
#' @param ind 
#' @param peptide 
#' @param tabdata 
#' @param metacond 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
meanImp_emmeans <- function(ind,peptide=1,tabdata,metacond) {
  tempdata <- cbind(reponse=tabdata[peptide,,ind],data.frame(fact=factor(metacond)))
  templm <- lm(reponse~fact,data=tempdata)
  return(summary(emmeans::emmeans(templm,specs = "fact"))$emmean)
}