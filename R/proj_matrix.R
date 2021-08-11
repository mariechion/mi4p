#' Title
#'
#' @param VarRubin.matrix 
#' @param metadata 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
proj_matrix <- function(VarRubin.matrix, metadata) {
  return(as.numeric(lapply(VarRubin.matrix, function(aaa){
    DesMat = DAPAR::make.design(metadata)
    max(diag(aaa)%*%t(DesMat)%*%DesMat)
  })))
}