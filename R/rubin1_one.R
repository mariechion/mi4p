# Computes first Rubin's rule for a given peptide

#' Title
#'
#' @param peptide 
#' @param data 
#' @param funcmean 
#' @param metacond 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
rubin1.one <- function(peptide,data,funcmean = meanImp_emmeans,metacond) {
  return(rowMeans(simplify2array(lapply(1:dim(data)[3],
                                        funcmean,
                                        tabdata=data,
                                        peptide=peptide,
                                        metacond=metacond))))
}