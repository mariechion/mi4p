norm.200.m100.sd1.vs.m200.sd1.list <- lapply(1:100, function(aaa){
  datasim <- data.frame(id.obs = 1:200, matrix(NA, 200, 10))
  datasim[1:10,-1] <- t(replicate(10,c(rnorm(5, 100, 1), rnorm(5, 200, 1))))
  datasim[-1:-10,-1] <- t(replicate(190,c(rnorm(5, 100, 1), rnorm(5, 100, 1))))
  return(datasim)
})

metadata <- data.frame(Sample.name = colnames(norm.200.m100.sd1.vs.m200.sd1.list[[1]][,-1]), 
                        Condition = as.factor(rep(c("A","B"), c(5,5))),
                        Bio.Rep = 1:ncol(norm.200.m100.sd1.vs.m200.sd1.list[[1]][,-1]))
