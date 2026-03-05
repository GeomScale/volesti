inner_ball_nloptr <- function(A, b){

  library(nloptr)

  d <- ncol(A)
  m <- nrow(A)

  # objective: maximize r  → minimize -r
  eval_f <- function(x){

    r <- x[d+1]

    return(-r)
  }

  # constraints: A_i^T x + ||A_i|| r <= b_i
  eval_g_ineq <- function(x){

    center <- x[1:d]
    r <- x[d+1]

    norms <- apply(A,1,function(v) sqrt(sum(v^2)))

    return(A %*% center + norms*r - b)
  }

  # initial guess
  x0 <- c(rep(0,d),0.1)

  res <- nloptr(
    x0=x0,
    eval_f=eval_f,
    eval_g_ineq=eval_g_ineq,
    opts=list(
      algorithm="NLOPT_LD_MMA",
      xtol_rel=1e-8
    )
  )

  center <- res$solution[1:d]
  radius <- res$solution[d+1]

  return(list(
    center=center,
    radius=radius,
    status=res$status
  ))
}