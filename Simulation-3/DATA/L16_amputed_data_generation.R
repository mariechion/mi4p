# --------------------------------------- #
#   Script for amputed data generation    #
# --------------------------------------- #

library("doMC")
registerDoMC(detectCores(all.tests = FALSE, logical = TRUE)-1)

# ---- RANDOM NUMBER GENERATION ---- #
set.seed(17)

# ---- LOAD FUNCTIONS ---- #
source("Functions/MV_generation.R")

# ---- LOAD COMPLETE DATA ---- #
load(file = "Simulation-3/DATA/L16_noMV_data")

# ---- MV 1% ---- #
L16.MV1pct.NA.data <- foreach(iforeach =  L16.noMV.data,
                          .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach, prop_NA = 0.01)

# ---- MV 5% ---- #
L16.MV5pct.NA.data <- foreach(iforeach =  L16.noMV.data,
                              .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach, prop_NA = 0.05)

# ---- MV 10% ---- #
L16.MV10pct.NA.data <- foreach(iforeach =  L16.noMV.data,
                              .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach, prop_NA = 0.10)

# ---- MV 15% ---- #
L16.MV15pct.NA.data <- foreach(iforeach =  L16.noMV.data,
                              .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach, prop_NA = 0.15)

# ---- MV 20% ---- #
L16.MV20pct.NA.data <- foreach(iforeach =  L16.noMV.data,
                              .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach, prop_NA = 0.20)

# ---- MV 25% ---- #
L16.MV25pct.NA.data <- foreach(iforeach =  L16.noMV.data,
                              .errorhandling = 'stop', .verbose = T) %dopar% 
  MVgen(dataset = iforeach, prop_NA = 0.25)
