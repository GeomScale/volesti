#' Apply a random rotation to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this function applies a random rotation.
#' 
#' @param P A convex polytope. It is an object from class (a) HPolytope or (b) VPolytope or (c) Zonotope.
#' 
#' @return A random rotation of the polytope that is given as an input. The return class is the same as the input class.
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' P = GenSimplex(2,'H')
#' listHpoly = rand_rotate(P)
#' 
#' # rotate a V-polytope (3d cube)
#' P = GenCube(3, 'V')
#' matVpoly = rand_rotate(P)
#' 
#' # rotate a 5-dimensional zonotope defined by the Minkowski sum of 15 segments
#' Z = GenZonotope(3,6)
#' MatZono = rand_rotate(Z)
#' @export
rand_rotate <- function(P){
  
  if (!missing(P)) {
    repr = class(P)[1]
    if (repr == "HPolytope") {
      vpoly = FALSE
      Zono = FALSE
    } else if(repr == "VPolytope") {
      vpoly = TRUE
      Zono = FALSE
    } else if(repr == "Zonotope") {
      vpoly = FALSE
      Zono = TRUE
    } else {
      stop("Not a known polytope representation.")
    }
    Mat = P$get_mat()
    dimension = dim(Mat)[2] - 1
    walk_length = 10 + floor( dimension / 10 )
  } else {
    stop("No polytope is given.")
  }
  
  Mat = rotating(Mat, Zono, vpoly)
  
  # get elements "matrix" and "vector"
  # remove first row
  Mat = Mat[-c(1),]
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  A = Mat[,-c(1)]
  if (vpoly) {
    PP = VPolytope(V=A)
  }else if (Zono) {
    PP = Zonotope(G=A)
  } else {
    PP = HPolytope("A"=A, "b"=b)
  }
  return(PP)
}
