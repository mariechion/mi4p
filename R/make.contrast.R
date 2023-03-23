

#' This function builds the contrast matrix
#' 
#' @title Builds the contrast matrix
#' 
#' @param design The data.frame which correspond to the pData function of 
#' MSnbase
#' 
#' @param condition xxxxx
#' 
#' @param contrast An integer that Indicates if the test consists of the 
#' comparison of each biological condition versus each of the other ones 
#' (Contrast=1; for example H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#' or each condition versus all others (Contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#'  H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#'  
#' @return A constrat matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek originally in 
#' the \code{DAPAR} package. Included in this package since \code{DAPAR} was to be removed from 
#' Bioconductor >= 3.15. 
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' design <- make.design(Biobase::pData(Exp1_R25_pept))
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' make.contrast(design, conds)
#' }
#' 
#' @export
#' 
make.contrast <- function(design, condition, contrast=1){
  
  
  #######################################################
  aggreg.column.design=function(design,Condition){
    nb.cond=length(levels(Condition))
    name.col=colnames(design)
    name.cond=NULL
    nb.col=NULL
    for (i in 1:nb.cond){
      col.select=NULL
      col.name.begin=paste("Condition",i, sep = "")
      nc=nchar(col.name.begin)
      for (j in 1:length(design[1,])){
        if (substr(name.col[j], 1, nc)==col.name.begin){
          col.select=c(col.select,j)
        }
      }
      name.aggreg=NULL
      for (j in 1:length(col.select)){
        name.aggreg=paste(name.aggreg,name.col[col.select[j]],sep="+")
      }
      name.aggreg=substr(name.aggreg, 2, nchar(name.aggreg))
      name.cond=c(name.cond,name.aggreg)
      nb.col=c(nb.col,length(col.select))
    }
    return(list(name.cond,nb.col))
  }
  
  
  
  nb.cond=length(levels(condition))
  r=aggreg.column.design(design,condition)
  label.agg=r[[1]]
  nb.agg=r[[2]]
  k=1
  
  if (contrast == 1){
    ## Contrast for One vs One
    contra=rep(0,sum(1:(nb.cond-1)))
    for (i in 1:(nb.cond-1)){
      for (j in (i+1):nb.cond){
        contra[k]=c(paste("(",label.agg[i],")/",
                          nb.agg[i],"-(",label.agg[j],")/",
                          nb.agg[j]))
        k=k+1
      }
    }
  } else if (contrast==2){
    ## Contrast for One vs All
    contra=rep(0,nb.cond)
    for (i in 1:(nb.cond)){
      contra[k]=c(paste("(",label.agg[i],")/",nb.agg[i]))
      nb=sum(nb.agg[(1:nb.cond)[(1:nb.cond)!=i]])
      for (j in (1:nb.cond)[(1:nb.cond)!=i]){
        contra[k]=c(paste(contra[k],"-(",label.agg[j],")/",nb))
      }
      k=k+1
    }
  }
  
  return(contra)
}
