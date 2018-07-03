#takes matrix .ine style and return matrix A and b: Ax<=b
modifyMat <- function(A){
  
  b=A[,1]
  A2=A[,-c(1)]
  retList=list("matrix"=A2, "vector"=b)
  return(retList)
  
}

#functiion to get a read.cs('path/to/file.ine') and return matrix for VolEsti()
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

#function to run the tests
testRvolEsti <- function(){
  path=getwd()
  path=paste0(substr(path, start=1, stop=nchar(path)-7),'/test/test_data/')
  print(path)
  listofexamples=list.files(path)
  
  for(i in 1:length(listofexamples)){
    x=read.csv(paste0(path,listofexamples[i]))
    print(listofexamples[i])
    A=ineToMatrix(x)
    VolEsti(list("matrix"=A,"verbose"=TRUE))
  }
  
}

#Take matrices A,b: Ax<=b and compute the chebychev center using lpSolveAPI library
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


#Main function. Takes a list and computes volume
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
  print(Cheb_ball)
  vol=vol_R(A,W,e,Cheb_ball,coordinate,rounding,verbose)
  #print(paste0('magnitude: ',ceiling(-log10(vol))))
  tim=proc.time()-tim
  if(verbose){
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(vol)
  
}
