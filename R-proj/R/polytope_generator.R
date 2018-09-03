
# This is an internal function. No Rd file.
polytope_generator <- function(Zono, Vpoly_gen, kind_gen, dim_gen, m_gen) {
  
  if (dim_gen<0) {
    print('Wrong Inputs.. You could read the documentation.')
    return(FALSE)
  }
  if (Zono && m_gen<0) {
    print('Wrong Inputs.. You could read the documentation.')
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
  sample_only = FALSE
  numpoints = 0
  variance = 0
  coordinate = TRUE
  rounding = FALSE
  verbose=FALSE
  #-------------------#
  
  Mat = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              Vpoly, Zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, 
              rotate_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  
  # get elements "matrix" and "vector"
  retList = modifyMat(Mat)
  if (Vpoly || Zono){
    # in V-polytope or Zonotope case return only the marix
    return(retList$matrix)
  } else {
    return(retList)
  }
  
}