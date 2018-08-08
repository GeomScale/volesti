#' Compute the Chebychev ball of a H-polytope.
#' 
#' For a H-polytope described by a \eqn{m\times d} matrix A and a d-dimensional vector b, s.t.: \eqn{Ax\leq b}, this function computes the largest inscribed ball of that polytope by solving the corresponding linear program.
#'
#' @param A the matrix of the H-polytope.
#' @param b The d-dimensional vector b that containes the constants of the facets.
#' @return A d+1-dimensional vector that containes the chebychev ball. The first d coordinates corresponds to the center and the last one to the radius of the chebychev ball.
#' @examples
#' #compute the Chebychev ball of a 2d unit simplex
#' A = matrix(c(-1,0,0,-1,1,1),ncol=2,nrow=3,byrow=TRUE)
#' b = c(0,0,1)
#' ball_vec = CheBall(A,b)
CheBall <- function(A,b){
  
  d=dim(A)[2]
  m=dim(A)[1]
  
  lprec <- make.lp(m, d+1)
  norm_row=rep(0,m)

  for(j in 1:m){
    norm_row[j]=norm(A[j,],type="2")
  }
  for(i in 1:d){
    set.column(lprec, i, A[,i])
  }
  set.column(lprec, d+1, norm_row)
  
  set.objfn(lprec, c(rep(0,d),c(-1)))
  set.constr.type(lprec, rep("<=",m))
  set.rhs(lprec, b)
  
  set.bounds(lprec, lower = rep(-Inf,d), columns = 1:d)
  
  solve(lprec)
  
  return(get.variables(lprec))
}

