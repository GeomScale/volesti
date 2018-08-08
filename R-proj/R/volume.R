#' The main R function for volume approximation of a convex H-Polytope
#'
#' @param list("argument"=value) A list that includes alla the parameters of the algorithm
#' @param path The path to an ine or ext file that describes the H-polytope. If path is given then "matrix" and "vector" inputs are not needed
#' @param matrix The matrix A of the polytope. If it is in ine format then the input "vector" is not needed
#' @param vector The vector b that containes the constants of the hyperplanes
#' @param Walk_length Optional. Declare the number of the steps for the random walk, default is 10+d/10
#' @param error Optional. Declare the goal for the approximation error. Default is 1 for volesti and 0.2 for CV.
#' @param Chebychev Optional. A d+1 vector that containes the chebychev center in the first d coordinates and the radius of the chebychev ball in the last coordinate
#' @param annealing Optional. A boolean parameter to use CV algorithm. Default value is false.
#' @param win_len Optional. The size of the window for the ratios' approximation in CV algorithm. Default value is win_len=4*(dimension^2)+500
#' @param C Optional. a constant for the upper boud of variance/mean^2 in schedule annealing
#' @param N optional. The number of points we sample in each step of schedule annealing in CV algorithm. Default value is N=500*C+(dimension^2)/2
#' @param ratio Optional. parameter of schedule annealing, larger ratio means larger steps in schedule annealing. Default value is ratio=1-1/dimension
#' @param frac Optional. the fraction of the total error to spend in the first gaussian. Default value is frac=0.1
#' @param ball_walk Optional. Boolean parameter to use ball walk, only for CV algorithm .Default value is False
#' @param delta Optional. The radius for the ball walk
#' @param verbose Optional. A boolean parameter for printing. Default is False
#' @param vpoly A boolean parameter, has to be true when a V-polytope is given as input
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is True
#' @param rounding Optional. A boolean parameter to activate the rounding option. Default value is False
#' @param test Optional. A boolean parameter. Declare if the current excecution is a test or not. Default value is False
#' @return The approximation of the volume of an H-polytope
#' @examples
#' VolEsti(list("path"=/path/to/ine/file, "verbose"=TRUE))
volume <- function(Inputs){
  
  Vpoly=FALSE
  if(!is.null(Inputs$Vpoly)){
    Vpoly = Inputs$Vpoly
  }
  if(!is.null(Inputs$path)){
    A=ineToMatrix(read.csv(Inputs$path))
    r=A[1,]
    A=A[-c(1),]
    x=modifyMat(A)
    A=x$matrix
    b=x$vector
  }else if(!is.null(Inputs$vector)){
    b=Inputs$vector
    A=-Inputs$matrix
    d=dim(A)[2]+1
    m=dim(A)[1]
    r=rep(0,d)
    r[1]=m
    r[2]=d
  }else if(!is.null(Inputs$matrix)){
    if(Vpoly){
      A=Inputs$matrix
      d=dim(A)[2]+1
      m=dim(A)[1]
      b=rep(1,m)
      r=rep(0,d)
      r[1]=m
      r[2]=d
    }else{
      r=Inputs$matrix[1,]
      Inputs$matrix=Inputs$matrix[-c(1),]
      x=modifyMat(Inputs$matrix)
      A=x$matrix
      b=x$vector
    }
  }else{
    if(Vpoly){
      print('No V-polytope defined from input!')
    }else{
      print('No H-polytope defined from input!')
    }
    return(-1)
  }
  A=matrix(cbind(b,A),ncol=dim(A)[2]+1)
  A=matrix(rbind(r,A),ncol=dim(A)[2])
  
  Cheb_ball=rep(0,dim(A)[2]+5)
  if(!is.null(Inputs$Chebychev)){
    Cheb_ball=Inputs$Chebychev
  }
  annealing=FALSE
  if(!is.null(Inputs$annealing)){
    annealing=Inputs$annealing
  }
  verbose=FALSE
  if(!is.null(Inputs$verbose)){
    verbose=Inputs$verbose
  }
  test=FALSE
  if(!is.null(Inputs$test)){
    test=Inputs$test
  }
  coordinate=TRUE
  if(!is.null(Inputs$coordinate)){
    coordinate=Inputs$coordinate
  }
  rounding=FALSE
  if(!is.null(Inputs$rounding)){
    rounding=Inputs$rounding
  }
  if(!is.null(Inputs$Walk_length)){
    W=Inputs$Walk_length
  }else{
    if(annealing){
      W=1
    }else{
      W=10+floor(dim(A)[2]/10)
    }
  }
  if(!is.null(Inputs$error)){
    e=Inputs$error
  }else{
    if(annealing){
      e=0.2
    }else{
      e=1
    }
  }
  dimension=dim(A)[2]-1
  win_len=4*(dimension^2)+500
  if(!is.null(Inputs$window_len)){
    win_len=Inputs$window_len
  }
  C=2
  if(!is.null(Inputs$C)){
    C=Inputs$C
  }
  ratio=1-1/dimension
  if(!is.null(Inputs$ratio)){
    ratio=Inputs$ratio
  }
  N=500*C+(dimension^2)/2
  if(!is.null(Inputs$N)){
    N=Inputs$N
  }
  frac=0.1
  if(!is.null(Inputs$frac)){
    frac=Inputs$frac
  }
  ball_walk=FALSE
  if(!is.null(Inputs$ball_walk)){
    ball_walk=Inputs$ball_walk
  }
  delta=-1
  if(!is.null(Inputs$delta)){
    delta=Inputs$delta
  }

  sample_only = FALSE
  variance = 0
  numpoints = 0
  tim=proc.time()
  vol=vol_R(A,W,e,Cheb_ball,annealing,win_len,N,C,ratio,frac,ball_walk,delta,Vpoly,sample_only,numpoints,variance,coordinate,rounding,verbose)
  tim=proc.time()-tim
  if(verbose || test){
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(vol[1,1])
  
}

