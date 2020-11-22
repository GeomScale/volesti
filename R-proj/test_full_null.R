library(volesti)
library(Matrix)
library(Rmosek)

path_mets='/home/tolis/data/metabolic_mat/iIT341.mat'

P = metabolic_net_2_polytope(path_mets)

pre_proc_list = fast_preprocess_with_mosek(P)

A = pre_proc_list$Aeq
m = dim(A)[1]
n=dim(A)[2]
N = null_space_SparseQR(pre_proc_list$row_ind, pre_proc_list$col_ind, pre_proc_list$values, m, n)

#print(N)

print(A%*%N)

print(dim(N))
