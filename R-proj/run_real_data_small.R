library(volesti)

#dm_covariance_matrices_70d <- readRDS("~/volume_approximation/R-proj/dm_covariance_matrices_70d.rds")
europe_ex_ch_covariance_matrices_small <- readRDS("~/temporal_repos/volume_approximation/R-proj/europe_ex_ch_covariance_matrices_small.rds")
N = length(europe_ex_ch_covariance_matrices_small$lCov)

all_samples_smallest_21_80 = matrix(list(), 0, 1)

n = dim(europe_ex_ch_covariance_matrices_small$lCov[[1]])[2] #number of assets

win = 40 #window length
M = 1000 #portfolios to generate for each level of volatility
something_went_wong_smallest = matrix(,0,2) #structure to store the instances where the method failed

#for (i in 125:(k-win+1)) {
for (i in 21:N) {
  
  vol_level_samples = matrix(list(), 0, 1)
  
  sigma = europe_ex_ch_covariance_matrices_small$lCov[[i]] #get the covariance matrix
  
  Cs = europe_ex_ch_covariance_matrices_small$lVola_targets[[i]] #get m level of volatilities
  m = length(Cs)
  
  for (j in 1:m) {
    c = Cs[j]
    print(paste0('i = ',as.character(i), ' j = ', as.character(j)))
    samples = sample_ptfs_constant_volatility(sigma, c, M) #sample M points, with c volatility when the cov. matrix is sigma
    
    correctness = check_correctness(samples, sigma, c) # check if the points lie in the simplex and have the requested value of volatility
    print(dim(samples))
    print(correctness)
    
    vol_level_samples[[length(vol_level_samples) + 1]] = samples
    
    if (!correctness) {
      something_went_wong_smallest = rbind(something_went_wong_smallest, c(i,j))
    }
  }
  all_samples_smallest_21_80[[length(all_samples_smallest_21_80) + 1]] = vol_level_samples
  saveRDS(all_samples_smallest_21_80, file = "all_samples_europe_smallest_21_80.rds")
  saveRDS(something_went_wong_smallest, file = "something_went_wong_smallest.rds")
}


