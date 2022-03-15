#' This function check the validity of the experimental design
#' 
#' @title Check if the design is valid
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A boolean
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek originally in 
#' the DAPAR package. Included in this package since DAPAR is to be removed from 
#' Bioconductor >= 3.15. 
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' check.design(Biobase::pData(Exp1_R25_pept)[,1:3])
#' }
#' 
#' @export
#' 
check.design <- function(sTab){
  res <- list(valid=FALSE,warn=NULL)
  
  names <- colnames(sTab)
  level.design <- ncol(sTab)-2
  
  
  res <- check.conditions(sTab$Condition)
  if (!res$valid){
    return(res)
  }
  # Check if all the column are fullfilled
  
  if (level.design == 1){
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)){
      res <- list(valid=FALSE,warn="The Bio.Rep colmumn are not full filled.")
      return(res)
    }
  }
  else if (level.design == 2){
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)){
      res <- list(valid=FALSE,warn="The Bio.Rep colmumn are not full filled.")
      return(res)
    }else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)){
      res <- list(valid=FALSE,warn="The Tech.Rep colmumn are not full filled.")
      return(res)
    }
  }
  else if (level.design == 3){
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)){
      res <- list(valid=FALSE,warn="The Bio.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)){
      res <- list(valid=FALSE,warn="The Tech.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Analyt.Rep) || (NA %in% sTab$Analyt.Rep)){
      res <- list(valid=FALSE,warn="The Analyt.Rep colmumn are not full filled.")
      return(res)
    }
  }
  
  # Check if the hierarchy of the design is correct
  if (level.design == 1){res <- test.design(sTab[,c("Condition", "Bio.Rep")])}
  else if (level.design == 2){res <- test.design(sTab[,c("Condition", "Bio.Rep","Tech.Rep")])}
  else if (level.design == 3){
    res <- test.design(sTab[,c("Condition", "Bio.Rep","Tech.Rep")])
    if (res$valid)
    {
      res <- test.design(sTab[,c("Bio.Rep","Tech.Rep", "Analyt.Rep")])
      
    }
  }
  
  return(res)
}

