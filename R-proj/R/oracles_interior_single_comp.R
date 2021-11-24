#' export
eval_f0 <- function(x, A, b) {
  return( (sum(x^2) - 1)^2 )
}


#' export
eval_g0 <- function(x, A, b) {
  q = A%*%x-b
  #if(sum(q>0)>0) {
  #  return(1)
  #} else {
  #  return(-1)
  #}
  return(q)
}


#' export
eval_grad_f0 <- function(x, A, b) {
  return( (sum(x^2) - 1)*x )
}


#' export
eval_jac_g0 <- function(x, A, b) {
  #d=length(x)
  return(A)
}


#' export
compute_interior_point_single_component <- function(A, b, x0) {
  
  bb = b - A%*%x0
  p = -x0
  d=length(x0)
  res0 <- nloptr::nloptr( x0=p,
                  eval_f=eval_f0,
                  eval_grad_f=eval_grad_f0,
                  lb = rep(-6, d),
                  ub = rep(6, d),
                  eval_g_ineq = eval_g0,
                  eval_jac_g_ineq = eval_jac_g0,
                  opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel"=1.0e-8, maxeval = 50000),
                  A = A,
                  b = bb )
  x = c(res0$solution) + c(x0)
  print(sum(eval_g0(x,A,b)>0))
  print(sqrt(sum((x-x0)^2)))
  #NLOPT_LD_SLSQP
  #NLOPT_LN_COBYLA
  return(x)
}
