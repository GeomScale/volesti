#' export
eval_f0_2 <- function(x, A, b) {
  return( -tail(x,1) )
}


#' export
eval_g0_2 <- function(x, A, b) {
  d = length(x)-1
  q = A%*%x-b
  w = sqrt(sum(x[1:d]^2))+x[d+1]-1
  q = c(q,w)
  #if(sum(q>0)>0) {
  #  return(1)
  #} else {
  #  return(-1)
  #}
  return(q)
}


#' export
eval_grad_f0_2 <- function(x, A, b) {
  d = length(x)
  grad_vec = rep(0,d)
  grad_vec[d] = -1
  return(grad_vec)
}


#' export
eval_jac_g0_2 <- function(x, A, b) {
  #d=length(x)
  return(A)
}


#' export
get_A <- function(A) {
  
  m = dim(A)[1]
  d = dim(A)[2]
  B = matrix(0,m,d+1)
  B[1:m,1:d] = A
  
  for (i in 1:m) {
    q = A[i,]
    norm_q = sqrt(sum(q^2))
    B[i,d+1] = norm_q
  }
  return(B)
}


#' export
compute_max_ball <- function(A, b, x0, x) {
  
  b = b - A%*%x0
  A = get_A(A)
  
  xx = x - x0
  p = c(c(xx), 0.00001)
  d=length(x0)
  res0 <- nloptr::nloptr( x0=p,
                          eval_f=eval_f0_2,
                          #eval_grad_f=eval_grad_f0_2,
                          #lb = rep(-600, d+1),
                          #ub = rep(600, d+1),
                          eval_g_ineq = eval_g0_2,
                          #eval_jac_g_ineq = eval_jac_g0_2,
                          opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel"=1.0e-5, maxeval = 10000),
                          A = A,
                          b = b )
  y = c(res0$solution)
  y[1:d] = y[1:d] + c(x0)
  #print(sum(eval_g0(x,A,b)>0))
  #print(sqrt(sum((x-x0)^2)))
  #NLOPT_LD_SLSQP
  #NLOPT_LN_COBYLA
  return(y)
}



#' export
get_point_on_component <- function(A, b, x0, V, Sind) {
  
  #n = length(Sind)
  d = length(x0)
  
  y = rep(0, d)
  
  #Xs = matrix(, d, 0)
  
  #for (i in 1:n) {
    
    ind_verts = Sind
    n_verts = length(ind_verts)
    rad_max = 0
    p=c()
    
    for (j in 1:n_verts) {
      
      v = V[,ind_verts[j]] - y
      res_int = ball_line_intersection(y, v, x0, 1)
      
      if (!res_int$intersect) {
        next
        #stop('no intersection')
      }
      if (is_in_component(y, x0, A, b, V, ind_verts, TRUE) && sqrt(sum((y-x0)^2)) > 1) {
        #print("isin")
        inds = 1:n
        inds = inds[-c(i)]
        ind_verts_temp = S_Vindices[[ inds[1] ]]
        v = V[, ind_verts_temp[1]] - y
        res_int = ball_line_intersection(y, v, x0, 1)
        #x = y + res_int$tmin[1]*v
        if (res_int$tmin[1] < 1 && res_int$tmin[1] > 0) {
          x = y + res_int$tmin[1]*v
        } else if(res_int$tmax[1] < 1 && res_int$tmax[1] > 0) {
          x = y + res_int$tmax[1]*v
        }
      } else if(sqrt(sum((y-x0)^2)) < 1) {
        #x = y + res_int$tmin[1]*v
        if (res_int$tmax[1] < 1 && res_int$tmax[1] > 0) {
          x = y + res_int$tmax[1]*v
        } else if(res_int$tmin[1] < 1 && res_int$tmin[1] > 0) {
          x = y + res_int$tminx[1]*v
        }
      } else {
        #x = y + res_int$tmin[1]*v
        if (res_int$tmax[1] < 1 && res_int$tmax[1] > 0) {
          x = y + res_int$tmax[1]*v
        } else if(res_int$tmin[1] < 1 && res_int$tmin[1] > 0) {
          x = y + res_int$tmin[1]*v
        }
      }
      
      if (!is_in_component(x, x0, A, b, V, ind_verts, TRUE)) {
        stop('point outside')
      }
      
      rad = get_rad_cap(x, A, b, x0)
      if(rad>rad_max) {
        rad_max = rad
        p = x
      }
    }
    
  
  return(p)
}


#' export
get_center_max_ball <- function(A, b, x0, V, Sind) {
  
  x = get_point_on_component(A, b, x0, V, Sind)
  v = x0 - x
  lambdas =  A %*% v / (b - A%*%x)
  l_max = max(lambdas)
  l_max = 1 / l_max
  x = x + (l_max/2)*v
  
  d = length(x)
  
  print("computing max ball..")
  y = compute_max_ball(A, b, x0, x)
  print("max ball computed..")
  xc = y[1:d]
  
  q = A%*%xc-b
  if(sum(q>0) > 0){
    stop("[max_ball] center outside simplex")
  }
  if (sqrt(sum((xc-x0)^2))+y[d+1]-1 > 0) {
    stop("[max_ball] too big interior ball")
  }
  
  return(xc)
}


#' export
get_points_on_component_imp <- function(A, b, x0, V, Sind) {
  
  xc = get_center_max_ball(A, b, x0, V, Sind)
  n = length(Sind)
  d = length(x0)
  
  
  
  Xs = matrix(, d, 0)
  
  #for (i in 1:n) {
    
    ind_verts = Sind
    n_verts = length(ind_verts)
    
    for (j in 1:n_verts) {
      y = rep(0, d)
      fail_num = 0
      rad_max = 0
      p=c()
      for (ii in 1:2) {
  
        v = V[,ind_verts[j]] - y
        res_int = ball_line_intersection(y, v, x0, 1)
      
        if (!res_int$intersect) {
          fail_num = fail_num + 1
          if (fail_num == 2) {
            stop('num_fails = 2')
          }
          y = xc
          next
        }
        if (is_in_component(y, x0, A, b, V, ind_verts, TRUE) && sqrt(sum((y-x0)^2)) > 1) {
          #print("isin")
          inds = 1:n
          inds = inds[-c(i)]
          ind_verts_temp = S_Vindices[[ inds[1] ]]
          v = V[, ind_verts_temp[1]] - y
          res_int = ball_line_intersection(y, v, x0, 1)
          #x = y + res_int$tmin[1]*v
          if (res_int$tmin[1] < 1 && res_int$tmin[1] > 0) {
            x = y + res_int$tmin[1]*v
          } else if(res_int$tmax[1] < 1 && res_int$tmax[1] > 0) {
            x = y + res_int$tmax[1]*v
          }
        } else if(sqrt(sum((y-x0)^2)) < 1) {
          #x = y + res_int$tmin[1]*v
          if (res_int$tmax[1] < 1 && res_int$tmax[1] > 0) {
            x = y + res_int$tmax[1]*v
          } else if(res_int$tmin[1] < 1 && res_int$tmin[1] > 0) {
            x = y + res_int$tminx[1]*v
          }
        } else {
          #x = y + res_int$tmin[1]*v
          if (res_int$tmax[1] < 1 && res_int$tmax[1] > 0) {
            x = y + res_int$tmax[1]*v
          } else if(res_int$tmin[1] < 1 && res_int$tmin[1] > 0) {
            x = y + res_int$tmin[1]*v
          }
        }
      
        if (!is_in_component(x, x0, A, b, V, ind_verts, TRUE)) {
          stop('point outside')
        }
      
        rad = get_rad_cap(x, A, b, x0)
        if(rad>rad_max) {
          rad_max = rad
          p = x
        }
        y = xc
      }
      Xs = cbind(Xs,p)
    }
    
  return(Xs)
}



