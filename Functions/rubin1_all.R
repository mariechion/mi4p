source("Functions/rubin1.one.R")

# Computes first Rubin's rule for all peptides
rubin1.all <- function(data, metacond, funcmean = meanImp_emmeans, is.parallel = F) {
  if (is.parallel) {
    res<-foreach(iforeach=1:dim(data)[1], .combine=cbind, 
                 .errorhandling = 'remove', .verbose = F) %dopar% 
      rubin1.one(iforeach,data=data,
                 funcmean=funcmean,metacond=metacond)
    res<-t(simplify2array(res))
    rownames(res) <- rownames(data)
  }
  else {
    res <- t(simplify2array(lapply(1:dim(data)[1],
                                   rubin1.one,
                                   data=data,
                                   funcmean=funcmean,
                                   metacond=metacond)))
  }
  return(res)
}
