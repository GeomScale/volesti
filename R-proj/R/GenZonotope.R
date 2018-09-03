

GenZonotope <- function(DimGen, NumGen) {
  
  Zono = TRUE
  kind_gen = 0
  Vpoly_gen = FALSE
  
  ListMat = polytope_generator(Zono, Vpoly_gen, kind_gen, DimGen, NumGen)
  
  return(ListMat)
  
}
  
  # set flag for verbose mode
  verbose=FALSE
  if(!is.null(Inputs$verbose)){
    verbose=Inputs$verbose
  }
  
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
  Vpoly_gen = FALSE
  #-------------------#
  
  Mat = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              Vpoly, Zono, gen_only, Vpoly_gen, kind_gen, DimGen, NumGen, round_only, 
              rotate_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  
  # get elements "matrix" and "vector"
  retList = modifyMat(Mat)
  return(retList$matrix)
}