#' Compute the Chebychev ball of a H-polytope.
#' 
#' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}, this function computes the largest inscribed ball (Chebychev ball) of that polytope by solving the corresponding linear program.
#' This function need suggested R-package lpSolveAPI.
#'
#' @param A The matrix of the H-polytope.
#' @param b The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets.
#' @return A \eqn{(d+1)}-dimensional vector that containes the Chebychev ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the Chebychev ball.
#' @examples
#' # compute the Chebychev ball of a 2d unit simplex
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' ball_vec = CheBall(A,b)
CheBall <- function(A,b){
  
  if (!requireNamespace("lpSolveAPI", quietly = TRUE)) {
    stop("Package \"lpSolveAPI\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  d = dim(A)[2]
  m = dim(A)[1]
  
  lprec <- make.lp(m, d+1)
  norm_row = rep(0,m)

  # get 2-norm of each row
  for(j in 1:m){
    norm_row[j] = norm(A[j,], type="2")
  }
  
  # set coeffs
  for (i in 1:d) {
    set.column(lprec, i, A[,i])
  }
  set.column(lprec, d+1, norm_row)
  
  # set objective function
  set.objfn(lprec, c(rep(0,d), c(-1)))
  
  # set the type of constraints
  set.constr.type(lprec, rep("<=",m))
  set.rhs(lprec, b)
  
  set.bounds(lprec, lower = rep(-Inf,d), columns = 1:d)
  
  solve(lprec)
  
  return(get.variables(lprec))
}
