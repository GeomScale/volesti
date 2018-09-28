

copula2 <- function(h, E, numSlices, N) {
  
  if(missing(h)) {
    stop('Only one family is given.')
  }
  
  if(missing(E)) {
    stop('Only one family is given.')
  }
  
  if (dim(E)[1]!=dim(E)[2]) {
    stop("Matrix E is not square.")
  }
  
  if (length(h)!=dim(E)[1]){
    stop("The two families are defined for different dimensions.")
  }
  
  if (length(h)<2 || dim(E)[1]<2) {
    stop("The dimension of the two families has to be greater than 2.")
  }
  
  slices = 100
  if(!missing(numSlices)) {
    slices = numSlices
    if(slices<=0) {
      stop("Number of slices has to be a positive integer.")
    }
  }
  n = 4000000
  if (!missing(N)) {
    n = N
    if(n<=0) {
      stop("Number of sampled points has to be a positive integer.")
    }
  }
  
  Mat = copula_hypEll(h, E, slices, n)
  return(Mat)
    
}
