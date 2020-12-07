library(volesti)

## give the path to your mat file that represents the model
path = '~/snap/volume_approximation/R-proj/metabolic_mat_files/e_coli_core.mat'

## request effectiveness n = 1000 and sample steady states
## If you would like to sample the Recon2D_v04 or the Recon3D_301 
## from https://www.vmh.life/ you have to set TRUE one of the corresponding flags 
result_list = generate_steady_states(path, n = 1000, Recon2D_v04 = FALSE, Recon3D_301 = FALSE)

## compute the PSRF of all marginals
psrfs = psrf_univariate(result_list$samples)

## approximate the flux distribution of the reaction "Acetate kinase"
h1 = hist(result_list$steady_states[12,], 
     main="Acetate kinase", 
     xlab="Flux (mmol/gDW/h)", 
     border="black", 
     col="red", 
     xlim=c(min(result_list$steady_states[12,]), max(result_list$steady_states[12,])), 
     las=1, 
     breaks=50, 
     prob = TRUE)


## we could use the polytope of the last phase to sample more steady states as follows

N = 5000
inner_ball = get_max_inner_ball(result_list$HP_rounded$A, result_list$HP_rounded$b)
more_samples = sample_points(result_list$HP_rounded, n = N, 
                             random_walk = list("walk" = "aBiW", "walk_length" = 1,
                                                "starting_point"=inner_ball$center, 
                                                "L" = 6*sqrt(result_list$HP_rounded$dimension)*inner_ball$radius))

## map the points back to P0
samples_in_P0 = result_list$T %*% more_samples + 
  kronecker(matrix(1, 1, N), matrix(result_list$T_shift, ncol = 1))

## compute the steady states
more_steady_states = result_list$N %*% samples_in_P0 + 
  kronecker(matrix(1, 1, N), matrix(result_list$N_shift, ncol = 1))

