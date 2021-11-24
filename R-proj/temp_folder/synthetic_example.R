library(volesti)


#n = dim(Ret)[2] #number of assets

n = 100 #for synthetic data

#k = dim(Ret)[1] #number of weekly returns
win = 40 #window length
M = 5000 #portfolios to generate for each level of volatility
m = 5 #levels of volatility for time period
parameters = get_parameters(n-1) #parameters of the simulated annealing algorithm
something_went_wong = matrix(,0,2) #structure to store the instances where the method failed


sigma = rWishart::rWishart(1, 100, diag(n), covariance = TRUE)[, , 1] #sample a covariance from wishart distribution

#we need a new function to generate the sequence of volatilities
Cs = get_sequence_of_volatilities_2(sigma, m) #get m level of volatilities

c=Cs[1]

nu = parameters$nu
lb = parameters$lb
ub = parameters$ub
Nu = parameters$Nu
W = parameters$W
psrf_target = parameters$psrf_target
error = parameters$error
WW = parameters$WW

n <- dim(sigma)[2]

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

x0=center_2

Xs = get_points_on_components(A, b, center_2, V, S, S_Vindices, V_ind_out)

X = sample_component_psrf(A,b,Xs[,1],10000,1,x0,1.2)

Sind = S_Vindices[[1]]
res = get_max_cap(X, A, b, x0, V, Sind, dim(X)[2])

mu = res$center
sigma = cov(t(X))


vol = compute_fischer_volume(A,b,mu,sigma,x0,X,10,0.1)




