GenVpoly <- function(dimension, NumGen) {
  
  kind_gen = -2
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, NumGen)
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = VPolytope(V = Mat, t=2)
  
  return(P)
  
}
