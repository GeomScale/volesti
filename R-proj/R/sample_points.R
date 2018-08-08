

sample_points <- function(Inputs){
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
  
  numpoints=1
  if(!is.null(Inputs$numpoints)){
    numpoints=Inputs$numpoints
  }
  
  Cheb_ball=rep(0,dim(A)[2]+5)
  if(!is.null(Inputs$Chebychev)){
    Cheb_ball=Inputs$Chebychev
  }
  gaussian=FALSE
  if(!is.null(Inputs$gaussian)){
    gaussian=Inputs$gaussian
  }
  variance=1
  if(!is.null(Inputs$variance)){
    variance=Inputs$variance
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
  W=10+floor((dim(A)[2]-1)/10)
  if(!is.null(Inputs$Walk_length)){
    W=Inputs$Walk_length
  }
  ball_walk=FALSE
  if(!is.null(Inputs$ball_walk)){
    ball_walk=Inputs$ball_walk
  }
  delta=-1
  if(!is.null(Inputs$delta)){
    delta=Inputs$delta
  }
  rounding=FALSE
  e=0
  win_len=0
  C=0
  ratio=0
  N=0
  frac=0
  
  
  sample_only = TRUE
  tim=proc.time()
  points=vol_R(A,W,e,Cheb_ball,gaussian,win_len,N,C,ratio,frac,ball_walk,delta,Vpoly,sample_only,numpoints,variance,coordinate,rounding,verbose)
  tim=proc.time()-tim
  if(verbose || test){
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(points)
  
}