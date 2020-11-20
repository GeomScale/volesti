null_space_and_shift <- function(row_ind, col_ind, values, Aeq, beq) {
  
  m = dim(Aeq)[1]
  n = dim(Aeq)[2]
  print('solving system')
  x = solve_undetermined_system_lu(Aeq, beq)
  #x = solve_overdetermined_linear_system(row_ind, col_ind, values, beq, m, n)
  print('system solved')
  
  NN = pracma::nullspace(Aeq)
  
  res_list = list()
  res_list$N_shift = x
  res_list$N = NN
  
  return(res_list)
}
