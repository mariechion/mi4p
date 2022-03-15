#' This function builds the design matrix 
#' 
#' @title Builds the design matrix
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' make.design(Biobase::pData(Exp1_R25_pept))
#' }
#' 
#' @export
make.design <- function(sTab){
  
  if (!check.design(sTab)$valid){
    warning("The design matrix is not correct.")
    warning(check.design(sTab)$warn)
    return(NULL)
  }
  
  n <- ncol(sTab)
  if (n==1 || n==2){
    stop("Error in design matrix dimensions which must have at least 3 columns.")
  }
  
  res <- do.call(paste0("make.design.", (n-2)),list(sTab))
  
  return(res)
}



#' This function builds the design matrix for design of level 1
#' 
#' @title Builds the design matrix for designs of level 1
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' make.design.1(Biobase::pData(Exp1_R25_pept))
#' }
#' 
#' @export
#' 
# make.design.1 <- function(sTab){
#
#   Conditions <- factor(sTab$Condition,  levels=unique(sTab$Condition))
#   nb_cond=length(unique(Conditions))
#   nb_samples <- nrow(sTab)
#
#   #CGet the number of replicates per condition
#   # nb_Rep_decal=rep(0,nb_cond_decal)
#   # for (i in 1:nb_cond_decal){
#   #   nb_Rep_decal[i]=sum((Conditions_decal==unique(Conditions_decal)[i]))
#   # }
#
#   design=matrix(0,nb_samples,nb_cond)
#   #n0_decal=1
#   coln=NULL
#   for (j in 1:nb_cond){
#     test <- rep
#     coln=c(coln,paste("Condition",j,collapse=NULL,sep=""))
#     design[,j]=as.integer(Conditions==unique(Conditions)[j])
#   }
#   colnames(design)=coln
#
#   return(design)
# }


make.design.1 <- function(sTab){
  Conditions <- factor(sTab$Condition, ordered = TRUE)
  nb_cond=length(unique(Conditions))
  nb_samples <- nrow(sTab)
  
  #CGet the number of replicates per condition
  nb_Rep=rep(0,nb_cond)
  for (i in 1:nb_cond){
    nb_Rep[i]=sum((Conditions==unique(Conditions)[i]))
  }
  
  design=matrix(0,nb_samples,nb_cond)
  n0=1
  coln=NULL
  for (j in 1:nb_cond){
    coln=c(coln,paste("Condition",j,collapse=NULL,sep=""))
    design[(n0:(n0+nb_Rep[j]-1)),j]=rep(1,length((n0:(n0+nb_Rep[j]-1))))
    n0=n0+nb_Rep[j]
  }
  colnames(design)=coln
  
  return(design)
  
}




#' This function builds the design matrix for design of level 2
#' 
#' @title Builds the design matrix for designs of level 2
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' make.design.2(Biobase::pData(Exp1_R25_pept))
#' }
#' 
#' @importFrom stats model.matrix rnorm
#' 
#' @export
#' 
make.design.2=function(sTab){
  Condition <- factor(sTab$Condition, levels=unique(sTab$Condition))
  RepBio <- factor(sTab$Bio.Rep, levels=unique(sTab$Bio.Rep))
  
  #Renome the levels of factor
  levels(Condition)=c(1:length(levels(Condition)))
  levels(RepBio)=c(1:length(levels(RepBio)))
  
  #Initial design matrix
  df <- rep(0,nrow(sTab))
  names(df) <- rownames(sTab)
  design=model.matrix(df~0+Condition:RepBio)
  
  #Remove empty columns in the design matrix
  design=design[,(apply(design,2,sum)>0)]
  #Remove identical columns in the design matrix
  coldel=-1
  for (i in 1:(length(design[1,])-1)){
    d2=as.matrix(design[,(i+1):length(design[1,])]);
    for (j in 1:length(d2[1,])){
      d2[,j]=d2[,j]-design[,i];
    }
    e=as.matrix(rnorm(length(design[,1]),10,1));
    sd2=t(e)%*%d2
    liste=which(sd2==0)
    coldel=c(coldel,liste+i)
  }
  design=design[,(1:length(design[1,]))!=coldel]
  colnames(design)=make.names(colnames(design))
  return(design)
}




#' This function builds the design matrix for design of level 3
#' 
#' @title Builds the design matrix for designs of level 3
#' 
#' @param sTab The data.frame which correspond to the pData function of MSnbase
#' 
#' @return A design matrix
#' 
#' @author Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek originally in 
#' the DAPAR package. Included in this package since DAPAR is to be removed from 
#' Bioconductor >= 3.15. 
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' sTab <-cbind(Biobase::pData(Exp1_R25_pept), Tech.Rep=1:6)
#' make.design.3(sTab)
#' }
#' 
#' @importFrom stats model.matrix rnorm
#' 
#' @export
#' 
make.design.3 <- function(sTab){
  
  Condition <- factor(sTab$Condition, levels=unique(sTab$Condition))
  RepBio <- factor(sTab$Bio.Rep, levels=unique(sTab$Bio.Rep))
  RepTech <- factor(sTab$Tech.Rep, levels=unique(sTab$Tech.Rep))
  
  
  #Rename the levels of factor
  levels(Condition)=c(1:length(levels(Condition)))
  levels(RepBio)=c(1:length(levels(RepBio)))
  levels(RepTech)=c(1:length(levels(RepTech)))
  
  
  #Initial design matrix
  df <- rep(0,nrow(sTab))
  names(df) <- rownames(sTab)
  design=model.matrix(df~0+Condition:RepBio:RepTech)
  
  #Remove empty columns in the design matrix
  design=design[,(apply(design,2,sum)>0)]
  
  #Remove identical columns in the design matrix
  coldel=-1
  for (i in 1:(length(design[1,])-1)){
    d2=as.matrix(design[,(i+1):length(design[1,])]);
    for (j in 1:length(d2[1,])){
      d2[,j]=d2[,j]-design[,i];
    }
    e=as.matrix(rnorm(length(design[,1]),10,1));
    sd2=t(e)%*%d2
    liste=which(sd2==0)
    coldel=c(coldel,liste+i)
  }
  design=design[,(1:length(design[1,]))!=coldel]
  colnames(design)=make.names(colnames(design))
  return(design)
}


