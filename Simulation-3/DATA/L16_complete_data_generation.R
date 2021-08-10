# ------------------------------------------------- #
#   Script for complete data generation             #
#   Simulated dataset after Lazar et al. (2016)     #
# ------------------------------------------------- #

set.seed(17)


# 1000 peptides
# 20% differentially abaundant i.e. 200
# 20 replicates, split into 2 groups of 10 replicates, for 2 conditions

# Model (Simplified version of Karpievitch et al. 2012)
# y_ij = P_i + G_ik + epsilon_ij

# "Pi is randomly generated from a Gaussian distribution with mean mu  = 1.5
# and standard deviation sigma = 0.5"
P_i <- rnorm(n = 1000, mean = 1.5, sd = 0.5)

# "Gik randomly drawn from the distribution previously mentioned"
G_ik <- c(rnorm(n = 200, mean = 1.5, sd = 0.5), rep(0,800))

# "the random error term has also been simulated by random draws from 
# a Gaussian distribution with zero mean and standard deviation sigmaϵ = 0.5."
eps_ij <- rnorm(n = 1000, mean = 0, sd = 0.5)

# Final dataset
cond1 <- replicate(n = 10, expr = rnorm(n = 1000, mean = 1.5, sd = 0.5) + 
                                   rnorm(n = 1000, mean = 0, sd = 0.5))
cond2 <- replicate(n = 10, expr = rnorm(n = 1000, mean = 1.5, sd = 0.5) + 
                                   c(rnorm(n = 200, mean = 1.5, sd = 0.5), rep(0,800)) +
                                   rnorm(n = 1000, mean = 0, sd = 0.5))
sim.data <- as.data.frame(cbind(cond1,cond2))

# Simulation of 100 datasets
L16.noMV.data <- lapply(1:100, function(aaa){
  cond1 <- replicate(n = 10, expr = rnorm(n = 1000, mean = 1.5, sd = 0.5) + 
                       rnorm(n = 1000, mean = 0, sd = 0.5))
  cond2 <- replicate(n = 10, expr = rnorm(n = 1000, mean = 1.5, sd = 0.5) + 
                       c(rnorm(n = 200, mean = 1.5, sd = 0.5), rep(0,800)) +
                       rnorm(n = 1000, mean = 0, sd = 0.5))
  sim.data <- as.data.frame(cbind(cond1,cond2))
  return(sim.data)
})

# Metadata
L16_metadata <- data.frame(Sample.name = as.factor(colnames(sim.data.list[[1]])),
                          Condition = as.factor(rep(c("A","B"),c(10,10))),
                          Bio.Rep = 1:20)

