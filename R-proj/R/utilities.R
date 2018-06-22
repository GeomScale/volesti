modifyMat <- function(A){
  
  b=A[,1]
  A2=A[,-c(1)]
  retList=list("matrix"=A2, "vector"=b)
  return(retList)
  
}


CheBall <- function(A,b){
  
  d=dim(A)[2]
  m=dim(A)[1]
  
  lprec <- make.lp(m, d+1)
  norm_row=rep(0,m)
  #A2=A
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
  #set.bounds(lprec, lower = 0, columns = d+1)
  #set.bounds(lprec, upper = rep(Inf,d+1), columns = 1:(d+1))
  
  solve(lprec)
  #ret=c(c(get.variables(lprec)),c(get.objective(lprec)))
  
  
  return(get.variables(lprec))
  
}


VolEsti <- function(Inputs){
  
  #A=Inputs$matrix
  if(!is.null(Inputs$vector)){
    b=Inputs$vector
    A=-Inputs$matrix
  }else{
    r=Inputs$matrix[1,]
    Inputs$matrix=Inputs$matrix[-c(1),]
    x=modifyMat(Inputs$matrix)
    A=x$matrix
    b=x$vector
  }
  
  if(!is.null(Inputs$Chebychev)){
    xc=Inputs$Chebychev
  }else{
    xc=CheBall(A,b)
  }
  #print(A)
  #print(b)
  A=matrix(cbind(b,A),ncol=dim(A)[2]+1)
  #print(A)
  A=matrix(rbind(r,A),ncol=dim(A)[2])
  #return(list("matrix"=A,"vector"=b,"cheb"=xc))
  return(vol_R(A,10,1,xc))
  
}
