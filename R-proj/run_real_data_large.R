library(volesti)

#dm_covariance_matrices_70d <- readRDS("~/volume_approximation/R-proj/dm_covariance_matrices_70d.rds")
usa_covariance_matrices_large <- readRDS("~/volume_approximation/R-proj/usa_covariance_matrices_large.rds")
N = length(usa_covariance_matrices_large$lCov)
print(N)

all_samples_largest_1_81 = matrix(list(), 0, 1)

n = dim(usa_covariance_matrices_large$lCov[[1]])[2] #number of assets

M = 5000 #portfolios to generate for each level of volatility
something_went_wong_smallest = matrix(,0,2) #structure to store the instances where the method failed
ignore_smallest_components = FALSE

start_index = 56
end_index = 76

for (index in seq(from=start_index, to=end_index, by=4)) {
  
  all_samples_largest_1_81 = matrix(list(), 0, 1)

  for (i in index:(index+9)) {
  
    vol_level_samples = matrix(list(), 0, 1)
  
    sigma = usa_covariance_matrices_large$lCov[[i]] #get the covariance matrix
  
    Cs = usa_covariance_matrices_large$lVola_targets[[i]] #get m level of volatilities
    m = length(Cs)
  
    for (j in 1:m) {
      c = Cs[j]
      print(paste0('i = ',as.character(i), ' j = ', as.character(j)))
      res = sample_ptfs_constant_volatility(sigma, c, M, ignore_smallest_components) #sample M points, with c volatility when the cov. matrix is sigma
    
      nn = length(res$overall_samples)
      print(paste0("nn = ",as.character(nn)))
      for (ii in 1:nn) {
        indxs = res$num_verts
        #print(indxs)
        samples = res$overall_samples[[ii]]
        correctness = check_correctness(samples, sigma, c) # check if the points lie in the simplex and have the requested value of volatility
        print(dim(samples))
        print(correctness)
      }
    
      vol_level_samples[[length(vol_level_samples) + 1]] = res
    
      if (!correctness) {
        something_went_wong_smallest = rbind(something_went_wong_smallest, c(i,j))
      }
    }
    all_samples_largest_1_81[[length(all_samples_largest_1_81) + 1]] = vol_level_samples
    saveRDS(all_samples_largest_1_81, file = paste0("all_samples_usa_largest_",as.character(index),"_",as.character(index+9),".rds"))
    saveRDS(something_went_wong_smallest, file = paste0("something_went_wong_usa_largest_",as.character(index),"_",as.character(index+9),".rds"))
  }

}
