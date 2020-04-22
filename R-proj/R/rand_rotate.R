#' Apply a random rotation to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this function applies a random rotation.
#' 
#' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
#' 
#' @return A list that contains the rotated polytope and the matrix of the linear transformation.
#'
#' @details Let \eqn{P} be the given polytope and \eqn{Q} the rotated one and \eqn{T} be the matrix of the linear transformation. 
#' \itemize{
#' \item{If \eqn{P} is in H-representation and \eqn{A} is the matrix that contains the normal vectors of the facets of \eqn{Q} then \eqn{AT} contains the normal vactors of the facets of \eqn{P}.}
#' \item{If \eqn{P} is in V-representation and \eqn{V} is the matrix that contains column-wise the vertices of \eqn{Q} then \eqn{T^TV} contains the vertices of \eqn{P}.}
#' \item{If \eqn{P} is a zonotope and \eqn{G} is the matrix that contains column-wise the generators of \eqn{Q} then \eqn{T^TG} contains the generators of \eqn{P}.}
#' \item{If \eqn{M} is a matrix that contains column-wise points in \eqn{Q} then \eqn{T^TM} contains points in \eqn{P}.}
#' }
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' P = GenSimplex(2,'H')
#' poly_matrix_list = rand_rotate(P)
#' 
#' # rotate a V-polytope (3d cube)
#' P = GenCube(3, 'V')
#' poly_matrix_list = rand_rotate(P)
#' 
#' # rotate a 5-dimensional zonotope defined by the Minkowski sum of 15 segments
#' Z = GenZonotope(3,6)
#' poly_matrix_list = rand_rotate(Z)
#' @export
rand_rotate <- function(P){
  
  #call rcpp rotating function
  Mat = rotating(P)
  
  n = P$dimension
  m=dim(Mat)[2]-n
  Tr = Mat[,-c(1:(dim(Mat)[2]-n))]
  Tr = Tr[1:n, 1:n]
  Mat = t(Mat[,1:m])
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  #A = Mat[,-c(1)]
  type = P$type
  if (type == 2) {
    PP = Vpolytope$new(A)
  }else if (type == 3) {
    PP = Zonotope$new(A)
  } else {
    PP = Hpolytope$new(A, b)
  }
  return(list("P" = PP, "T" = Tr))
}
