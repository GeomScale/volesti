GenHpoly <- function(dimension, NumGen) {
  
  kind_gen = -1
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, NumGen)
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = HPolytope(A = Mat, b = b, t=1)
  
  return(P)
  
}
