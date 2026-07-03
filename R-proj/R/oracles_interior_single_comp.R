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
  d=length(x0)
  p = rep(0,d) - x0
  
  rad = Inf
  
  m = nrow(A)
  for (i in 1:m) {
    q = A[i,]
    rad_temp = abs(b[i]) / sqrt(sum(q^2))
    if (rad_temp < rad) {
      rad = rad_temp
    }
  }
  #print(paste0('rad = ',as.character(rad)))
    
  while(TRUE) {
    
    p = boundary_randsphere(d,1)[,1]
    U = runif(1, 0, 1)
    U = U^(1/d)
    p = (rad*U)*c(p) - c(x0)
    
    y = A%*%p - bb
    if (sum(y>0) > 0) {
      stop('sampled point outside')
    }
    
    lb_coord = -10
    ub_coord = 10
    if (max(p) > ub_coord){
      ub_coord = 2*max(p)
    }
    if (min(p) < lb_coord){
      lb_coord = 2*min(p)
    }
  
    #print(max(p))
    #print(min(p))
  
  
    res0 <- nloptr::nloptr( x0=p,
                  eval_f=eval_f0,
                  eval_grad_f=eval_grad_f0,
                  lb = rep(lb_coord, d),
                  ub = rep(ub_coord, d),
                  eval_g_ineq = eval_g0,
                  eval_jac_g_ineq = eval_jac_g0,
                  opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel"=1.0e-8, maxeval = 50000),
                  A = A,
                  b = bb )
    x = c(res0$solution) + c(x0)
    #print(sum(eval_g0(x,A,b)>0))
    num_facets_valid = sum(eval_g0(x,A,b)<0)
    if (num_facets_valid != dim(A)[1]) {
      next
      print('interior point outside [single component case!')
    } else {
      break
    }
  }
  #stop('stop')
  #print(sum(eval_g0(x,A,b)==0))
  #rad = sqrt(sum((x-x0)^2))
  #NLOPT_LD_SLSQP
  #NLOPT_LN_COBYLA
  return(x)
}



###export
#compute_interior_point_single_component <- function(A, b, x0) {
  
#  x = compute_interior_point_single_component_preprocess(A, b, x0)
  
#}

