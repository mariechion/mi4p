## ----setup, include = FALSE---------------------------------------------------
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=FALSE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5
)

## ---- eval = FALSE------------------------------------------------------------
#  devtools::install_github("mariechion/mi4p")

## -----------------------------------------------------------------------------
library(mi4p)

## -----------------------------------------------------------------------------
set.seed(4619)
datasim <- protdatasim()
str(datasim)

## -----------------------------------------------------------------------------
attr(datasim, "metadata")

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
MV1pct.NA.data <- MVgen(dataset = datasim[,-1], prop_NA = 0.01)
MV1pct.NA.data

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
MV1pct.impMLE <- multi.impute(data = MV1pct.NA.data, conditions = attr(datasim,"metadata")$Condition, method = "MLE", parallel = FALSE)

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
print(paste(Sys.time(), "Dataset", 1, "out of", 1))
MV1pct.impMLE.VarRubin.Mat <- rubin2.all(data = MV1pct.impMLE, metacond = attr(datasim, "metadata")$Condition) 

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
print(paste("Dataset", 1, "out of",1, Sys.time()))
MV1pct.impMLE.VarRubin.S2 <- as.numeric(lapply(MV1pct.impMLE.VarRubin.Mat, function(aaa){
    DesMat = DAPAR::make.design(attr(datasim, "metadata"))
    return(max(diag(aaa)%*%t(DesMat)%*%DesMat))
  }))

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
MV1pct.impMLE.mi4limma.res <- mi4limma(qData = apply(MV1pct.impMLE,1:2,mean), 
                 sTab = attr(datasim, "metadata"), 
                 VarRubin = sqrt(MV1pct.impMLE.VarRubin.S2))
MV1pct.impMLE.mi4limma.res
(simplify2array(MV1pct.impMLE.mi4limma.res)$P_Value.A_vs_B_pval)[1:10]

(simplify2array(MV1pct.impMLE.mi4limma.res)$P_Value.A_vs_B_pval)[11:200]<=0.05

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
sum((simplify2array(MV1pct.impMLE.mi4limma.res)$P_Value.A_vs_B_pval)[1:10]<=0.05)/10

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
sum((simplify2array(MV1pct.impMLE.mi4limma.res)$P_Value.A_vs_B_pval)[11:200]<=0.05)/190

## ---- cache=TRUE, eval=LOCAL--------------------------------------------------
MV1pct.impMLE.dapar.res <-limmaCompleteTest.mod(qData = apply(MV1pct.impMLE,1:2,mean), sTab = attr(datasim, "metadata"))
MV1pct.impMLE.dapar.res

## -----------------------------------------------------------------------------
set.seed(4619)
norm.200.m100.sd1.vs.m200.sd1.list <- lapply(1:100, protdatasim)
metadata <- attr(norm.200.m100.sd1.vs.m200.sd1.list[[1]],"metadata")

## ---- eval=FALSE--------------------------------------------------------------
#  library(foreach)
#  doParallel::registerDoParallel(cores=NULL)
#  requireNamespace("foreach",quietly = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  MV1pct.NA.data <- foreach::foreach(iforeach =  norm.200.m100.sd1.vs.m200.sd1.list,
#                            .errorhandling = 'stop', .verbose = T) %dopar%
#    MVgen(dataset = iforeach[,-1], prop_NA = 0.01)

## ---- eval=FALSE--------------------------------------------------------------
#  MV1pct.impMLE <- foreach::foreach(iforeach =  MV1pct.NA.data,
#                           .errorhandling = 'stop', .verbose = F) %dopar%
#    multi.impute(data = iforeach, conditions = metadata$Condition,
#                 method = "MLE", parallel  = F)

## ---- eval=FALSE--------------------------------------------------------------
#  MV1pct.impMLE.VarRubin.Mat <- lapply(1:length(MV1pct.impMLE), function(index){
#    print(paste(Sys.time(), "Dataset", index, "out of", length(MV1pct.impMLE)))
#    rubin2.all(data = MV1pct.impMLE[[index]], metacond = metadata$Condition)
#  })

## ---- eval=FALSE--------------------------------------------------------------
#  MV1pct.impMLE.VarRubin.S2 <- lapply(1:length(MV1pct.impMLE.VarRubin.Mat), function(id.dataset){
#    print(paste("Dataset", id.dataset, "out of",length(MV1pct.impMLE.VarRubin.Mat), Sys.time()))
#    as.numeric(lapply(MV1pct.impMLE.VarRubin.Mat[[id.dataset]], function(aaa){
#      DesMat = DAPAR::make.design(metadata)
#      return(max(diag(aaa)%*%t(DesMat)%*%DesMat))
#    }))
#  })

## ---- eval=FALSE--------------------------------------------------------------
#  MV1pct.impMLE.mi4limma.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = T) %dopar%
#    mi4limma(qData = apply(MV1pct.impMLE[[iforeach]],1:2,mean),
#                   sTab = metadata,
#                   VarRubin = sqrt(MV1pct.impMLE.VarRubin.S2[[iforeach]]))
#  
#  MV1pct.impMLE.dapar.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = T) %dopar%
#    limmaCompleteTest.mod(qData = apply(MV1pct.impMLE[[iforeach]],1:2,mean),
#                          sTab = metadata)

