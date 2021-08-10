source("Functions/meanImp_emmeans.R")

# Computes first Rubin's rule for a given peptide
rubin1.one <- function(peptide,data,funcmean = meanImp_emmeans,metacond) {
  return(rowMeans(simplify2array(lapply(1:dim(data)[3],
                                        funcmean,
                                        tabdata=data,
                                        peptide=peptide,
                                        metacond=metacond))))
}