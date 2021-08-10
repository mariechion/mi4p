proj_matrix <- function(VarRubin.matrix, metadata) {
  return(as.numeric(lapply(VarRubin.matrix, function(aaa){
    DesMat = DAPAR::make.design(metadata)
    max(diag(aaa)%*%t(DesMat)%*%DesMat)
  })))
}