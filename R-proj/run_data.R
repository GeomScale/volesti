library(volesti)

msci_ci <- readRDS("~/volume_approximation/R-proj/data/msci_ci.rds")
Ret = msci_ci$X_est

n = dim(Ret)[2] #number of assets

#n = 100 #for synthetic data

k = dim(Ret)[1] #number of weekly returns
win = 40 #window length
M = 5000 #portfolios to generate for each level of volatility
m = 5 #levels of volatility for time period
parameters = get_parameters(n-1) #parameters of the simulated annealing algorithm
something_went_wong = matrix(,0,2) #structure to store the instances where the method failed

for (i in 1:(k-win+1)) {
#for (i in 1:1) {
  
  R = Ret[i:(i+win-1),] #consider the returns of this sliding window
  sigma = cov.qis(R) #compute the covariance matrix
  
  #sigma = rWishart::rWishart(1, 100, diag(n), covariance = TRUE)[, , 1] #sample a covariance from wishart distribution
  
  #we need a new function to generate the sequence of volatilities
  Cs = get_sequence_of_volatilities_2(sigma, m) #get m level of volatilities
  
  for (j in 1:m) {
    c = Cs[j]
    print(paste0('i = ',as.character(i), ' j = ', as.character(j)))
    samples = sample_ptfs_constant_volatility(sigma, c, M, parameters) #sample M points, with c volatility when the cov. matrix is sigma
    correctness = check_correctness(samples, sigma, c) # check if the points lie in the simplex and have the requested value of volatility
    print(dim(samples))
    print(correctness)
    if (!correctness) {
      something_went_wong = rbind(something_went_wong, c(i,j))
    }
  }
  
}

