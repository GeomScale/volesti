#' Compute the Chebychev ball of a H-polytope, P:= Ax<=b
#'
#' @param A the matrix of the H-polytope
#' @param b the vector with the constants of the hyperplanes
#' @return The Chebychev center of the Polytope discribed by the matrix \code{A} and the vector \code{b}
#' @examples
#' CheBall(A,b)
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

