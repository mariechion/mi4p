source("Functions/rubin1_one.R")

# Computes the between-imputation component in the 2nd Rubin's rule
# for a given peptide

rubin2bt.one <- function(peptide,data,funcmean,metacond){
  funcmean_p = rubin1.one(peptide=peptide,data=data,funcmean = funcmean, metacond=metacond)
  outer_mat_prod = function(ind,peptide,funcmean_peptide){return(matrix((funcmean(ind,peptide=peptide,tabdata=data, metacond=metacond)-funcmean_peptide),ncol=1)%*%(funcmean(ind,peptide=peptide,tabdata=data, metacond=metacond)-funcmean_peptide))}
  return(1/(dim(data)[3]-1)*apply(simplify2array(lapply(1:dim(data)[3],outer_mat_prod,peptide=peptide,funcmean_p)),c(1,2),sum))    
}