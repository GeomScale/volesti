null_space_and_shift <- function(Aeq, beq) {
  
  x = solve_undetermined_system_lu(Aeq, beq)
  
  NN = pracma::nullspace(Aeq)
  
  res_list = list()
  res_list$N_shift = x
  res_list$N = NN
  
  return(res_list)
}
