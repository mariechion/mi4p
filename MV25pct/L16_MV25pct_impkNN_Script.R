library("doMC")
registerDoMC(detectCores(all.tests = FALSE, logical = TRUE)-1)

# Load functions
source("Functions/multi_impute.R")
source("Functions/rubin2_all.R")
source("Functions/proj_matrix.R")
source("Functions/mi4limma.R")
source("Functions/limmaCompleteTest_mod.R")


# Charger la liste des 100 jeux de données amputés
load(file = "Simulation-3/DATA/L16_MV25pct_NA_data")
# Charger les métadonnées
load(file = "Simulation-3/DATA/L16_metadata")

# ---- RANDOM NUMBER GENERATION ---- #
set.seed(17)

# ---- IMPUTATION ---- #
L16.MV25pct.impkNN <- foreach(iforeach =  L16.MV25pct.NA.data,
                         .errorhandling = 'stop', .verbose = F) %dopar% 
  multi.impute(data = iforeach, conditions = L16_metadata$Condition, 
               parallel  = F, method = "kNN")

# ---- ESTIMATION ---- #
L16.MV25pct.impkNN.VarRubin.Mat <- lapply(1:length(L16.MV25pct.impkNN), function(index){
  print(paste(Sys.time(), "Dataset", index, "out of", length(L16.MV25pct.impkNN)))
  rubin2.all(data = L16.MV25pct.impkNN[[index]], metacond = L16_metadata$Condition) 
})

# ---- PROJECTION ---- #
L16.MV25pct.impkNN.VarRubin.S2 <- lapply(1:length(L16.MV25pct.impkNN.VarRubin.Mat), function(id.dataset){
  print(paste("Dataset", id.dataset, "out of",length(L16.MV25pct.impkNN.VarRubin.Mat), Sys.time()))
  as.numeric(lapply(L16.MV25pct.impkNN.VarRubin.Mat[[id.dataset]], function(aaa){
    DesMat = DAPAR::make.design(L16_metadata)
    return(max(diag(aaa)%*%t(DesMat)%*%DesMat))
  }))
}) 

# ---- MODERATED T-TEST ---- #
L16.MV25pct.impkNN.mi4limma.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = F) %dopar%
  mi4limma(qData = apply(L16.MV25pct.impkNN[[iforeach]],1:2,mean), 
                 sTab = L16_metadata, 
                 VarRubin = sqrt(L16.MV25pct.impkNN.VarRubin.S2[[iforeach]]))

L16.MV25pct.impkNN.dapar.res <- foreach(iforeach =  1:100,  .errorhandling = 'stop', .verbose = F) %dopar%
  limmaCompleteTest.mod(qData = apply(L16.MV25pct.impkNN[[iforeach]],1:2,mean), 
                        sTab = L16_metadata)