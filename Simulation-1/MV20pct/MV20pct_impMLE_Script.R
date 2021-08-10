library("doMC")
registerDoMC(detectCores(all.tests = FALSE, logical = TRUE)-1)

# Load functions
source("Functions/multi_impute.R")
source("Functions/rubin2_all.R")
source("Functions/proj_matrix.R")
source("Functions/mi4limma.R")
source("Functions/limmaCompleteTest_mod.R")

# Load simulated datasets
load("Simulation-1/DATA/norm_200_m100_sd1_vs_m200_sd1_list.RData")
# Load metadata
load("Simulation-1/DATA/norm_200_m100_sd1_vs_m200_sd1_metadata")

# ---- RANDOM NUMBER GENERATION ---- #
set.seed(17)

# ---- AMPUTATION ---- #
MV20pct.NA.data <- foreach(iforeach =  norm.200.m100.sd1.vs.m200.sd1.list,
                          .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach[,-1], prop_NA = 0.2)

# ---- IMPUTATION ---- #
MV20pct.impMLE <- foreach(iforeach =  MV20pct.NA.data,
                         .errorhandling = 'stop', .verbose = F) %dopar% 
  multi.impute(data = iforeach, conditions = metadata$Condition, 
               method = "MLE", parallel  = F)

# ---- ESTIMATION ---- #
MV20pct.impMLE.VarRubin.Mat <- lapply(1:length(MV20pct.impMLE), function(index){
  print(paste(Sys.time(), "Dataset", index, "out of", length(MV20pct.impMLE)))
  rubin2.all(data = MV20pct.impMLE[[index]], metacond = metadata$Condition) 
})

# ---- PROJECTION ---- #
MV20pct.impMLE.VarRubin.S2 <- lapply(1:length(MV20pct.impMLE.VarRubin.Mat), function(id.dataset){
  print(paste("Dataset", id.dataset, "out of",length(MV20pct.impMLE.VarRubin.Mat), Sys.time()))
  as.numeric(lapply(MV20pct.impMLE.VarRubin.Mat[[id.dataset]], function(aaa){
    DesMat = DAPAR::make.design(metadata)
    return(max(diag(aaa)%*%t(DesMat)%*%DesMat))
  }))
}) 

# ---- MODERATED T-TEST ---- #
MV20pct.impMLE.mi4limma.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = T) %dopar%
  mi4limma(qData = apply(MV20pct.impMLE[[iforeach]],1:2,mean), 
                 sTab = metadata, 
                 VarRubin = sqrt(MV20pct.impMLE.VarRubin.S2[[iforeach]]))

MV20pct.impMLE.dapar.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = T) %dopar%
  limmaCompleteTest.mod(qData = apply(MV20pct.impMLE[[iforeach]],1:2,mean), 
                        sTab = metadata)
