library(volesti)

# build the 100-dimensional coanonical simplex
d = 100
A = -diag(d)
b = rep(0,d)
Aeq = matrix(1, 1, d)
beq = c(1)

## request effective sample size ess = 1000 and sample
result_list = samples_uniform_portfolios(A=A, b=b, Aeq=Aeq, beq=beq, ess = 1000)

## compute the PSRF of all marginals in the full dimensional space
psrfs = psrf_univariate(result_list$samples)

## approximate the marginal distribution of the 12-th marginal (weight)
h1 = hist(result_list$random_portfolios[12,], 
     main="Asset name", 
     xlab="Weight", 
     border="black", 
     col="red", 
     xlim=c(min(result_list$random_portfolios[12,]), max(result_list$random_portfolios[12,])), 
     las=1, 
     breaks=50, 
     prob = TRUE)


## use the polytope of the last phase to sample more portfolios
N = 5000
more_random_portfolios = sample_from_last_polytope(result_list, N)

## to compute a better density estimation
## approximate the marginal distribution
## using the total number of portfolios that we have generated
h2 = hist(c(result_list$random_portfolios[12,], more_random_portfolios[12, ]), 
          main="Asset name", 
          xlab="Weight", 
          border="black", 
          col="blue", 
          xlim=c(min(c(result_list$random_portfolios[12,], more_random_portfolios[12, ])), max(c(result_list$random_portfolios[12,], more_random_portfolios[12, ]))), 
          las=1, 
          breaks=50, 
          prob = TRUE)
