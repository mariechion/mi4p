library("doMC")
registerDoMC(detectCores(all.tests = FALSE, logical = TRUE)-1)

# Load functions
source("Functions/multi_impute.R")
source("Functions/rubin2_all.R")
source("Functions/proj_matrix.R")
source("Functions/mi4limma.R")
source("Functions/limmaCompleteTest_mod.R")


# Charger la liste des 100 jeux de données amputés
load(file = "Simulation-2/DATA/L16bis_MV15pct_NA_data")
# Charger les métadonnées
load(file = "Simulation-2/DATA/L16_metadata")

# ---- RANDOM NUMBER GENERATION ---- #
set.seed(17)

# ---- IMPUTATION ---- #
L16bis.MV15pct.impMLE <- foreach(iforeach =  L16bis.MV15pct.NA.data,
                         .errorhandling = 'stop', .verbose = F) %dopar% 
  multi.impute(data = iforeach, conditions = L16_metadata$Condition, 
               parallel  = F, method = "MLE")

# ---- ESTIMATION ---- #
L16bis.MV15pct.impMLE.VarRubin.Mat <- lapply(1:length(L16bis.MV15pct.impMLE), function(index){
  print(paste(Sys.time(), "Dataset", index, "out of", length(L16bis.MV15pct.impMLE)))
  rubin2.all(data = L16bis.MV15pct.impMLE[[index]], metacond = L16_metadata$Condition) 
})

# ---- PROJECTION ---- #
L16bis.MV15pct.impMLE.VarRubin.S2 <- lapply(1:length(L16bis.MV15pct.impMLE.VarRubin.Mat), function(id.dataset){
  print(paste("Dataset", id.dataset, "out of",length(L16bis.MV15pct.impMLE.VarRubin.Mat), Sys.time()))
  as.numeric(lapply(L16bis.MV15pct.impMLE.VarRubin.Mat[[id.dataset]], function(aaa){
    DesMat = DAPAR::make.design(L16_metadata)
    return(max(diag(aaa)%*%t(DesMat)%*%DesMat))
  }))
}) 

# ---- MODERATED T-TEST ---- #
L16bis.MV15pct.impMLE.mi4limma.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = F) %dopar%
  mi4limma(qData = apply(L16bis.MV15pct.impMLE[[iforeach]],1:2,mean), 
                 sTab = L16_metadata, 
                 VarRubin = sqrt(L16bis.MV15pct.impMLE.VarRubin.S2[[iforeach]]))

L16bis.MV15pct.impMLE.dapar.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = F) %dopar%
  limmaCompleteTest.mod(qData = apply(L16bis.MV15pct.impMLE[[iforeach]],1:2,mean), 
                        sTab = L16_metadata)