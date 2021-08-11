# Computes the within-variance component in the 2nd Rubin's rule 
# for a given peptide

#' Title
#'
#' @param peptide 
#' @param data 
#' @param funcvar 
#' @param metacond 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
rubin2wt.one <- function(peptide,data,funcvar,metacond){
  return(apply(simplify2array(lapply(1:dim(data)[3],
                                     funcvar,peptide=peptide,
                                     data=data,metacond=metacond)),c(1,2),mean))
}