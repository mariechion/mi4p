#' @title Check if the design is valid
#' 
#' @description 
#' This function checks the validity of the conditions.
#' 
#' This function was included from the \code{\link[DAPAR]{check.conditions}} 
#' function in the DAPAR package, since DAPAR is to be removed from 
#' Bioconductor >= 3.15. 
#' 
#' @param conds A vector
#' 
#' @return A list
#' 
#' @author Samuel Wieczorek as the author of 
#' \code{\link[DAPAR]{check.conditions}}.
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' check.conditions(Biobase::pData(Exp1_R25_pept)$Condition)
#' }
#' 
#' @export
#' 
check.conditions <- function(conds){
  res <- list(valid=TRUE,warn=NULL)
  
  if (("" %in% conds) || (NA %in% conds)){
    res <- list(valid=FALSE,warn="The conditions are note full filled.")
    return(res)
  }
  
  # Check if there is at least two conditions
  if (length(unique(conds)) < 2){
    res <- list(valid=FALSE,warn="The design must contain at least two conditions.")
    return(res)
  }
  
  
  # check if each condition has at least two values
  nValPerCond <- unlist(lapply(unique(conds), function(x){length(conds[which(conds==x)])}))
  if (all(nValPerCond < 2)){
    res <- list(valid=FALSE,warn="The design must contain at least two values per condition.")
    return(res)
  }
  
  return(res)
}
