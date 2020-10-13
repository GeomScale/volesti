#' Apply a random rotation to a convex polytope (H-polytope, V-polytope, zonotope or intersection of two V-polytopes)
#' 
#' Given a convex H- or V- polytope or a zonotope or an intersection of two V-polytopes as input, this function applies (a) a random rotation or (b) a given rotation by an input matrix \eqn{T}.
#' 
#' @param P A convex polytope. It is an object from class (a) Hpolytope, (b) Vpolytope, (c) Zonotope, (d) intersection of two V-polytopes.
#' @param rotation A list that contains (a) the rotation matrix T and (b) the 'seed' to set a spesific seed for the number generator.
#' 
#' @return A list that contains the rotated polytope and the matrix \eqn{T} of the linear transformation.
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
#' P = gen_simplex(2, 'H')
#' poly_matrix_list = rotate_polytope(P)
#' 
#' # rotate a V-polytope (3d cube)
#' P = gen_cube(3, 'V')
#' poly_matrix_list = rotate_polytope(P)
#' 
#' # rotate a 5-dimensional zonotope defined by the Minkowski sum of 15 segments
#' Z = gen_rand_zonotope(3, 6)
#' poly_matrix_list = rotate_polytope(Z)
#' @export
rotate_polytope <- function(P, rotation = list()) {
  
  seed = NULL
  if (!is.null(rotation$seed)) {
    seed = rotation$seed
  }
  
  if (is.null(rotation$T)) {
    T = NULL
  } else {
    T = rotation$T
  }
  
  #call rcpp rotating function
  Mat = rotating(P, T, seed)
  
  type = P@type
  
  if (type == 'Vpolytope') {
    n = dim(P@V)[2]
  }else if (type == 'Zonotope') {
    n = dim(P@G)[2]
  } else {
    n = dim(P@A)[2]
  }
  
  m = dim(Mat)[2] - n
  Tr = Mat[, -c(1:(dim(Mat)[2]-n)), drop = FALSE]
  Tr = Tr[1:n, 1:n]
  Mat = t(Mat[, 1:m])
  
  # first column is the vector b
  b = Mat[, 1]
  
  # remove first column
  A = Mat[, -c(1), drop = FALSE]
  
  if (type == 'Vpolytope') {
    PP = Vpolytope(V = A)
  }else if (type == 'Zonotope') {
    PP = Zonotope(G = A)
  } else {
    PP = Hpolytope(A = A, b = b)
  }
  return(list("P" = PP, "T" = Tr))
}
