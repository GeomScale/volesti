library(volesti)

#dm_covariance_matrices_70d <- readRDS("~/volume_approximation/R-proj/dm_covariance_matrices_70d.rds")
europe_ex_ch_covariance_matrices_small <- readRDS("~/volume_approximation/R-proj/europe_ex_ch_covariance_matrices_small.rds")
N = length(europe_ex_ch_covariance_matrices_small$lCov)

all_samples_smallest_21_80 = matrix(list(), 0, 1)

n = dim(europe_ex_ch_covariance_matrices_small$lCov[[1]])[2] #number of assets

win = 40 #window length
M = 1000 #portfolios to generate for each level of volatility
something_went_wong_smallest = matrix(,0,2) #structure to store the instances where the method failed

i=1
j=5

sigma = europe_ex_ch_covariance_matrices_small$lCov[[i]] #get the covariance matrix

Cs = europe_ex_ch_covariance_matrices_small$lVola_targets[[i]] #get m level of volatilities

c = Cs[j]

n <- dim(sigma)[2]

parameters = get_parameters(n-1)

nu = parameters$nu
lb = parameters$lb
ub = parameters$ub
Nu = parameters$Nu
W = parameters$W
psrf_target = parameters$psrf_target
error = parameters$error
WW = parameters$WW

A <- -diag(n)
b <- rep(0, 1)

Aeq <- matrix(rep(1, n), nrow = 1, ncol = n)
beq <- c(1)

N = pracma::nullspace(Aeq)
x0 <- rep(1, n)/n

b = b - A %*% x0
A = A %*% N
V = diag(n)
V = V - kronecker(matrix(1, 1, n), matrix(x0, ncol = 1))
V = t(N) %*% V
x0 = rep(0, n-1)

sigma_proj = t(N) %*% sigma %*% N
center = -MASS::ginv(sigma_proj) %*% (t(N) %*% sigma) %*% (rep(1,n)/n)

R = c + t(center)%*%sigma_proj%*%center - (rep(1,n)/n)%*%sigma%*%(rep(1,n)/n)
R = R[1]
sigma_proj = sigma_proj / R

b = b - A %*% center # shift ellipsoid's center to origin
x0 = x0 - center
V = V - kronecker(matrix(1, 1, n), matrix(center, ncol = 1))
center_2 = rep(0, n-1)

Tinv = t(chol(Matrix::nearPD(MASS::ginv(sigma_proj))$mat))
T = MASS::ginv(Tinv)

A = A %*% Tinv
x0_2 = T %*% x0
V = T %*% V

b = b - A %*% x0_2
center_2 = center_2 - x0_2
V = V - kronecker(matrix(1, 1, n), matrix(x0_2, ncol = 1))

components_res=find_components_new_vertices(V, center_2)

S = components_res[[1]]
S_Vindices = components_res[[2]]

if(length(components_res) > 2){
  V_ind_out = components_res[[3]]
} else {
  V_ind_out = c()
}

cmax = get_c_upper_bound(A, b, center_2)
cmin = 1

print(S_Vindices)

Xs = get_points_on_components_imp(A, b, center_2, V, S, S_Vindices)
xx = Xs[,1]
sample_component(A,b,xx,V,1000,1,center_2)

