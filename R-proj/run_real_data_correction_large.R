library(volesti)

#dm_covariance_matrices_70d <- readRDS("~/volume_approximation/R-proj/dm_covariance_matrices_70d.rds")
europe_ex_ch_covariance_matrices_large <- readRDS("~/volume_approximation/R-proj/europe_ex_ch_covariance_matrices_large.rds")
N = length(europe_ex_ch_covariance_matrices_large$lCov)

all_samples_largest_1_81 = matrix(list(), 0, 1)

n = dim(europe_ex_ch_covariance_matrices_large$lCov[[1]])[2] #number of assets

win = 40 #window length
M = 5000 #portfolios to generate for each level of volatility
something_went_wong_smallest = matrix(,0,2) #structure to store the instances where the method failed
with_samples = matrix(0, N, 5)


for (i in 18:N) {
  
  
  sigma = europe_ex_ch_covariance_matrices_large$lCov[[i]] #get the covariance matrix
  
  Cs = europe_ex_ch_covariance_matrices_large$lVola_targets[[i]] #get m level of volatilities
  m = length(Cs)
  
  for (j in 1:m) {
    c = Cs[j]
    print(paste0('i = ',as.character(i), ' j = ', as.character(j)))
    res = sample_ptfs_constant_volatility_to_correct(sigma, c, M) #sample M points, with c volatility when the cov. matrix is sigma
    
    if(res$computations) {
      with_samples[i,j]=1
      
      #correctness = check_correctness(res$samples, sigma, c) # check if the points lie in the simplex and have the requested value of volatility
      print(dim(res$samples))
      #print(correctness)
      #samples = res$samples
      saveRDS(res$samples, file = paste0("samples_europe_largest_corrected_",as.character(i),"_",as.character(j),".rds"))
      saveRDS(with_samples, file = "with_samples.rds")
      #res$samples = c()
      res = c()
      
    } 
    
  }
  
}

