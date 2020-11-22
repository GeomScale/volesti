library(volesti)
library(Matrix)
library(Rmosek)

path_mets='/home/tolis/data/metabolic_mat/iJO1366.mat'

P = metabolic_net_2_polytope(path_mets)
load('iJO1366_preprocess.RData')

print("Compute the null space to constraint")
rr = full_dimensional_polytope_with_arma(pre_proc_list$Aeq, pre_proc_list$beq)

print("Get full dimensional polytope")
A = P$A 
b = P$b
b = b - A %*% rr$N_shift
A = A %*% rr$N

m = dim(A)[1]
d = dim(A)[2]
rows_to_del = c()
for (i in 1:m) {
  if (sqrt(sum(A[i,]^2)) < 1e-06) {
    rows_to_del = c(rows_to_del, i)
  }
}
if (length(rows_to_del) > 0) {
  A = A[-rows_to_del, ]
  b = b[-rows_to_del]
}
print(paste0(length(rows_to_del), " facets removed"))

print("Compute scaling for numerical stability")
sc = central_scaling(A, b)
A = sc$A
b = sc$b
T_scale = sc$T_scale
scale_shift = sc$scale_shift


print("Rounding the polytope")
T_total = diag(d)
T_shift_total = rep(0, d)

if (d > 400 & FALSE) {
  z=get_max_inner_ball(A, b)
  res_list = rounding_max_ellipsoid_step(A, b, z$center, z$radius)
  A = res_list$A
  b = res_list$b
  T_total = res_list$T
  T_shift_total = res_list$shift
}

HP = Hpolytope$new(A = A, b = b)
ret_list = rounding_isotropic(HP)
if (FALSE) {
  save(ret_list, file = paste0(path,"rounding_matrices.RData"))
}
T_total = T_total %*% ret_list$T
T_shift_total = T_shift_total + ret_list$T_shift

HP = Hpolytope$new(A = ret_list$A, b = ret_list$b)
z=get_max_inner_ball(HP$A, HP$b)

N = 3000
print("Sample points from full diemnsional polytope")
samples =  sample_points(HP, random_walk = list("walk" = "aBiW", "starting_point" = z$center,
                                                "walk_length" = 1, "L" = 4*sqrt(d)*z$radius), n = N)

samples = T_total %*% samples + 
  kronecker(matrix(1, 1, N), matrix(T_shift_total, ncol = 1))

samples = T_scale %*% samples + 
  kronecker(matrix(1, 1, N), matrix(scale_shift, ncol = 1))

steady_states = rr$N %*% samples + 
  kronecker(matrix(1, 1, N), matrix(rr$N_shift, ncol = 1))

result_list = list()
result_list$HP_rounded = HP
result_list$steady_states = steady_states
result_list$N = rr$N
result_list$N_shift = rr$N_shift
result_list$T = ret_list$T
result_list$T_shift = ret_list$T_shift
result_list$minFluxes = pre_proc_list$minFluxes
result_list$maxFluxes = pre_proc_list$maxFluxes

save(result_list, file = "iJO1366_results.RData")

