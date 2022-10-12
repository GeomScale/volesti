library(volesti)

## Get Ax<b and Aeq x = beq from a metabolic model just to illustrate an example 
path = 'metabolic_mat_files/e_coli_core.mat'

P = metabolic_net_2_polytope(path, FALSE, FALSE)

## request effective sample size ess = 1000 and sample
result_list = samples_uniform_portfolios(A=P$A, b=P$b, Aeq=P$Aeq, beq=P$beq, ess = 1000)

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


## we could use the polytope of the last phase to sample more portfolios as follows
N = 5000
inner_ball = get_max_inner_ball(result_list$HP_rounded$A, result_list$HP_rounded$b)
more_samples = sample_points(result_list$HP_rounded, n = N, 
                             random_walk = list("walk" = "aBiW", "walk_length" = 1,
                                                "starting_point"=inner_ball$center, 
                                                "L" = 6*sqrt(result_list$HP_rounded$dimension)*inner_ball$radius))

## map the points back to P0
samples_in_P0 = result_list$T %*% more_samples + 
  kronecker(matrix(1, 1, N), matrix(result_list$T_shift, ncol = 1))

## compute the portfolios
more_random_portfolios = result_list$N %*% samples_in_P0 + 
  kronecker(matrix(1, 1, N), matrix(result_list$N_shift, ncol = 1))


## to compute a better density estimation
## approximate the marginal distribution
## using the total number of portfolios that we have generated
h2 = hist(c(result_list$more_random_portfolios[12,], more_random_portfolios[12, ]), 
          main="Asset name", 
          xlab="Weight", 
          border="black", 
          col="blue", 
          xlim=c(min(c(result_list$random_portfolios[12,], more_random_portfolios[12, ])), max(c(result_list$random_portfolios[12,], more_random_portfolios[12, ]))), 
          las=1, 
          breaks=50, 
          prob = TRUE)
