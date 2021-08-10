
source("Functions/rubin2bt_one.R")

# Computes the between-imputation component in the 2nd Rubin's rule
# for all peptides

rubin2bt.all <- function(data,funcmean = meanImp_emmeans, metacond, is.parallel = F){
  if (is.parallel){
    res <- foreach(iforeach=1:dim(data)[1], .errorhandling = 'remove', .verbose = verbose) %dopar% 
      rubin2bt.one(iforeach,data=data,funcmean=funcmean,metacond=metacond)
  }
  else {
    res <- lapply(1:dim(data)[1],rubin2bt.one,data=data,funcmean=funcmean,metacond=metacond)
  }
  return(res)
}
