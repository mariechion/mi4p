source("Functions/within_variance_comp_emmeans.R")

# Computes the within-variance component in the 2nd Rubin's rule 
# for a given peptide

rubin2wt.one <- function(peptide,data,funcvar,metacond){
  return(apply(simplify2array(lapply(1:dim(data)[3],
                                     funcvar,peptide=peptide,
                                     data=data,metacond=metacond)),c(1,2),mean))
}