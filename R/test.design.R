#' This function check xxxxx
#' 
#' @title Check if xxxxxx
#' 
#' @param tab A data.frame which correspond to xxxxxx
#' 
#' @return A list of two items
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek originally in 
#' the \code{DAPAR} package. Included in this package since \code{DAPAR} was to be removed from 
#' Bioconductor >= 3.15. 
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' test.design(Biobase::pData(Exp1_R25_pept)[,1:3])
#' }
#' 
#' @export
#' 
test.design <- function(tab){
  valid <- TRUE
  txt <- NULL
  level <- NULL
  
  level.a  <- factor(tab[,1], ordered = TRUE)
  level.b <- factor(tab[,2], ordered = TRUE)
  name.level.a <- colnames(tab)[1]
  name.level.b <-colnames(tab)[2]
  
  level.c <- NULL
  level.c <- if(ncol(tab)==3) factor(tab[,3], ordered = TRUE) else NULL
  name.level.c <- if(ncol(tab)==3) colnames(tab)[3] else NULL
  
  
  # verification intersection sur B
  ##verification de la non redondance'intersection vide entre les groupes
  uniqueA <- unique(level.a)
  ll <- lapply(uniqueA, function(x){as.character(level.b)[which(level.a==x)]})
  n <- NULL
  for (i in 1:(length(uniqueA)-1)){
    for (j in (i+1):length(uniqueA)){
      n <- c(n,intersect(ll[[i]], ll[[j]]))
    }
  }
  if (length(n) > 0){
    valid <- FALSE
    txt <- c(txt,paste0("The value ", n, " in column '", colnames(tab)[2], "' is not correctly set.\n"))
  }
  
  
  #verification si niveau hierarchique inf
  if (length(levels(level.a)) == length(levels(level.b))){
    ## c'est un design de niveau n-1 en fait
    valid <- FALSE
    txt <- c(txt,paste0("The column ",name.level.b, " is not informative. Thus, the design is not of level (n-1).\n"))
  }
  else if (!is.null(level.c)){
    if (length(levels(level.b)) == length(levels(level.c))){
      ## c'est un design de niveau n-1 en fait
      valid <- FALSE
      txt <- c(txt,paste0("The column ",name.level.c, " is not informative. Thus, the design is of level (n-1).\n"))
    }
  }
  
  #verification si niveau non informatif
  return(list(valid=valid, warn=txt))
}
