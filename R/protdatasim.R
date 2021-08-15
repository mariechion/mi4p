#' @title Data simulation function
#' 
#' @description Function to simulate benchmark datasets.
#'
#' @param iii A parameter useful to loop over for simulated lists of datasets. It has no effect.
#' @param nobs Number of peptides
#' @param nobs1 Number of peptides with differential expressions between the two conditions
#' @param ng1 Number of biological replicates in condition A
#' @param ng2 Number of biological replicates in condition B
#' @param mg1 Mean in condition A
#' @param mg2 Mean in condition B
#' @param dispg1 Dispersion in condition A
#' @param dispg2 Dispersion in condition B
#'
#' @return A data frame with the simulated and attribute metadata.
#' @export
#'
#' @examples
#' data_sim <- protdatasim()
#' attr(data_sim,"metadata")
#' 
#' norm.200.m100.sd1.vs.m200.sd1_list <- lapply(1:100, protdatasim)
#' attr(norm.200.m100.sd1.vs.m200.sd1_list[[1]],"metadata")
#' 
protdatasim <- function(iii=1,nobs=200,nobs1=10,ng1=5,ng2=5,mg1=100,mg2=200,dispg1=1,dispg2=1){
  datasim <- data.frame(id.obs = 1:nobs, matrix(NA, nobs, ng1+ng2))
  datasim[1:10,-1] <- t(replicate(nobs1,c(rnorm(ng1, mg1, dispg1), rnorm(ng2, mg2, dispg2))))
  datasim[-1:-10,-1] <- t(replicate(nobs-nobs1,c(rnorm(ng1, mg1, dispg1), rnorm(ng2, mg1, dispg1))))
  attr(datasim,"metadata") <- data.frame(Sample.name = colnames(datasim[,-1]), 
                                         Condition = as.factor(rep(c("A","B"), c(ng1,ng2))),
                                         Bio.Rep = 1:ncol(datasim[,-1]))
  return(datasim)
}
