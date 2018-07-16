#' takes a numerical matrix in ine format and return numerical matrix A and vector b: Ax<=b
#'
#' @param A the numerical matrix in ine format of the H-polytope
#' @return numerical matrix A and vector b: Ax<=b
#'
#' @examples
#' modifyMat(A)
modifyMat <- function(A){
  
  b=A[,1]
  A2=A[,-c(1)]
  retList=list("matrix"=A2, "vector"=b)
  return(retList)
  
}


#' functiion to get a ine file and return matrix A in ine format for VolEsti()
#'
#' @param read.cs('path/to/file.ine') the ine file desrcibing the H-polytope
#' @return The numerical matrix in ine format of \code{read.cs('path/to/file.ine')}
#' @examples
#' ineToMatrix(read.cs('path/to/data/cube40.ine'))
ineToMatrix <- function(P){
  
  r=as.character(P[3,1])
  count_sp=1
  str=""
  beg=0
  for(j in 1:nchar(r)){
    if(substr(r, start=j, stop=j)==" "){
      beg=beg+1
    }else{
      break
    }
  }
  for(i in seq(from=beg+1, to=nchar(r), by=1)){
    if(substr(r, start=i, stop=i)==" "){
      if(count_sp==1){
        m=as.numeric(str)
        str=""
        count_sp=count_sp+1
      }else{
        d=as.numeric(str)
        str=""
        break
      }
    }else{
      str=paste0(str,substr(r, start=i, stop=i))
    }
  }
  A=rep(0,d)
  A[1]=m
  A[2]=d
  newrow=rep(0,d)
  for(i in 4:(dim(P)[1]-2)){
    r=P[i,1]
    r=as.character(r)
    str=""
    count=1
    #print(r)
    beg=0
    for(j in 1:nchar(r)){
      if(substr(r, start=j, stop=j)==" "){
        beg=beg+1
      }else{
        break
      }
    }
    sp_bef=FALSE
    for(j in seq(from=beg+1, to=nchar(r), by=1)){
      
      if (substr(r, start=j, stop=j)==" "){
        if(sp_bef){
          next
        }
        sp_bef=TRUE
        newrow[count]=as.numeric(str)
        str=""
        count=count+1
      }else{
        str=paste0(str,substr(r, start=j, stop=j))
        sp_bef=FALSE
      }
    }
    A=rbind(A,newrow)
    newrow=rep(0,d)
  }
  A=matrix(A,ncol=dim(A)[2])
  
  return(A)
}

#' Run some experiments
#'
#' @param empty No inputs
#' @return Print the computed volumes and the total time
#' @examples
#' testRvolEsti()
testRvolEsti <- function(){
  path=getwd()
  path=paste0(substr(path, start=1, stop=nchar(path)-7),'/data/')
  print(path)
  listofexamples=list.files(path)
  
  for(i in 1:length(listofexamples)){
    x=read.csv(paste0(path,listofexamples[i]))
    print(listofexamples[i])
    A=ineToMatrix(x)
    VolEsti(list("matrix"=A,"test"=TRUE))
  }
  
}

#' Compute the Chebychev ball of a H-polytope, P:= Ax<=b
#'
#' @param A the matrix of the H-polytope
#' @param b the vector with the constants of the hyperplanes
#' @return The Chebychev center of the Polytope discribed by the matrix \code{A} and the vector \code{b}
#' @examples
#' CheBall(A,b)
CheBall <- function(A,b){
  
  d=dim(A)[2]
  m=dim(A)[1]
  
  lprec <- make.lp(m, d+1)
  norm_row=rep(0,m)

  for(j in 1:m){
    norm_row[j]=norm(A[j,],type="2")
  }
  for(i in 1:d){
    set.column(lprec, i, A[,i])
  }
  set.column(lprec, d+1, norm_row)
  
  set.objfn(lprec, c(rep(0,d),c(-1)))
  set.constr.type(lprec, rep("<=",m))
  set.rhs(lprec, b)
  
  set.bounds(lprec, lower = rep(-Inf,d), columns = 1:d)
  
  solve(lprec)
  
  return(get.variables(lprec))
}


#' The main R function for volume approximation of a convex H-Polytope
#'
#' @param list("path","matrix","vector","Chebychev","verbose","coordinate","rounding","Walk_length","error","test")
#' @param path The path to the ine file that describes the H-polytope. If path is given then "matrix" and "vector" inputs are not needed
#' @param matrix The matrix A of the polytope. If it is in ine format then the input "vector" is not needed
#' @param vector The vector b that containes the constants of the hyperplanes
#' @param Chebychev Optional. A d+1 vector that containes the chebychev center in the first d coordinates and the radius of the chebychev ball in the last coordinate
#' @param verbose Optional. A boolean parameter for printing. Default is False
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default is True
#' @param rounding Optional. A boolean parameter to activate the rounding option. Default is False
#' @param Walk_length Optional. Declare the number of the steps for the random walk, default is 10+d/10
#' @param error Optional. Declare the goal for the approximation error. Default is 1
#' @param test Optional. A boolean parameter. Declare if the current excecution is a test or not. Default is False
#' @return The approximation of the volume of an H-polytope
#' @examples
#' VolEsti(list("path"=/path/to/ine/file, "verbose"=TRUE))
VolEsti <- function(Inputs){
  
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
  }else if(!is.null(Inputs$matrix)){
    r=Inputs$matrix[1,]
    Inputs$matrix=Inputs$matrix[-c(1),]
    x=modifyMat(Inputs$matrix)
    A=x$matrix
    b=x$vector
  }else{
    print('No H-polytope defined from input!')
    return(-1)
  }
  
  if(!is.null(Inputs$Chebychev)){
    Cheb_ball=Inputs$Chebychev
  }else{
    #Cheb_ball=CheBall(A,b)
    Cheb_ball=rep(0,dim(A)[2]+5)
  }
  verbose=FALSE
  if(!is.null(Inputs$verbose)){
    if(Inputs$verbose){
      verbose=TRUE
    }else{
      verbose=FALSE
    }
  }
  test=FALSE
  if(!is.null(Inputs$test)){
    if(Inputs$test){
      test=TRUE
    }else{
      test=FALSE
    }
  }
  coordinate=TRUE
  if(!is.null(Inputs$coordinate)){
    if(Inputs$coordinate){
      coordinate=TRUE
    }else{
      coordinate=FALSE
    }
  }
  rounding=FALSE
  if(!is.null(Inputs$rounding)){
    if(Inputs$rounding){
      rounding=TRUE
    }else{
      rounding=FALSE
    }
  }
  if(!is.null(Inputs$Walk_length)){
    W=Inputs$Walk_length
  }else{
    W=10+floor(dim(A)[2]/10)
  }
  if(!is.null(Inputs$error)){
    e=Inputs$error
  }else{
    e=1
  }

  A=matrix(cbind(b,A),ncol=dim(A)[2]+1)
  A=matrix(rbind(r,A),ncol=dim(A)[2])
  tim=proc.time()
  vol=vol_R(A,W,e,Cheb_ball,coordinate,rounding,verbose)
  #print(paste0('magnitude: ',ceiling(-log10(vol))))
  tim=proc.time()-tim
  if(verbose || test){
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(vol)
  
}
