#' Internal function to generate polytopes
#' 
#' This function is used by polytope generator functions. It is an internal function and it is not suggested to use it.
#' 
#' @param Zono A boolean parameter to declare if the generated polytope has to be zonotope or not.
#' @param repr A string parameter to declare the representation of the polytope. Use 'H' for H-representation, 'V' for V-representation and 'zontope' for zonotopes.
#' @param kind_gen An integer to declare the kind of the polytope. Use '0' for zonotopes, '1' for cubes, '2' for cross polytopes, '3' for simplices, '4' for product of two simplices and '5' for skinny cubes. See polytope generator functions for more details.
#' @param dim_gen An integer to declare the dimension of the polytope.
#' @param m_gen Only for zonotopes. An integer to declare the number of segments.
#' 
#' @return For H-polytopes the return value is a list that containes a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b s.t.: \eqn{Ax\leq b}. For V-polytopes and zonotopes the return value is a \eqn{m\times d} matrix that containes row-wise the \eqn{d}-dimensional vertices or segments respectively.
#' 
#' @examples 
#' # create a 5-dimensional zonotope that is defined by the Minkowski sum of 10 segments
#' ZonoMat = polytope_generator(TRUE, 'zonotope', 0, 5, 10)
#' 
#' # create a 20-dimensional unit simplex in V-representation
#' PolyMat = polytope_generator(FALSE, 'V', 1, 20, -1) 
polytope_generator <- function(Zono, repr, kind_gen, dim_gen, m_gen) {
  
  if (repr == "V") {
    Vpoly_gen = TRUE
  } else if (repr == "H" || repr == "zonotope") {
    Vpoly_gen = FALSE
  } else {
    print('The second argument has to be V for V-polytopes or H for H-polytopes')
    return(FALSE)
  }
  
  if (dim_gen<0) {
    print('Wrong Inputs.. The dimension can not be negative.')
    return(FALSE)
  }
  if (Zono && m_gen<0) {
    print('Wrong Inputs.. The number of zonotope generators can not be negative.')
    return(FALSE)
  }
  
  if (Vpoly_gen) {
    if (kind_gen == 4) {
      print('No product of simplices can be generated in V-representation.. You could read the documentation.')
      return(FALSE)
    }
    if (kind_gen == 5) {
      print('No skinny cube can be generated in V-representation.. You could read the documentation.')
      return(FALSE)
    }
  }
  
  gen_only = TRUE
  
  #---------------------#
  A = matrix(c(0,0))
  round_only = FALSE
  rotate_only = FALSE
  W = 0
  e = 0
  Cheb_ball = c(0)
  annealing = FALSE
  win_len = 0
  N = 0
  C = 0
  ratio = 0
  frac = 0
  ball_walk = FALSE
  delta =0
  Vpoly = FALSE
  exact_zono = FALSE
  sample_only = FALSE
  numpoints = 0
  variance = 0
  coordinate = TRUE
  rounding = FALSE
  verbose=FALSE
  ball_only = FALSE
  #-------------------#
  
  Mat = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              Vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, 
              rotate_only, ball_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  
  # get elements "matrix" and "vector"
  # remove first row
  Mat = Mat[-c(1),]
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  if (Vpoly_gen || Zono){
    # in V-polytope or Zonotope case return only the marix
    return(A)
  } else {
    retList = list("A"=A, "b"=b)
    return(retList)
  }
  
}