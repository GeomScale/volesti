library(volesti)

rm(list = ls())
load("~/volume_approximation/R-proj/d_workspace_2.RData")
#d=dim(A)[2]
#vol = compute_fischer_volume(A,b,mu,x0,10,0.7)
#res1 = compute_fischer_annealing(A,b,mu,x0,10)
#res2 = compute_fischer_ratios(A,b,mu,x0,res1$variances, res1$ratios, 1000+floor((d*d)/2),10,0.5)
a_min = tail(res1$variances,1)
a_max = tail(res2$variances,1)
ratio_min = prod(1/res1$ratios)
ratio_max = prod(1/res2$ratios)
q = estimate_related_ratio_fischer(A, b, Xs[, 1], center_2, a_min, a_max, ratio_min, ratio_max, 10, 0.7)
