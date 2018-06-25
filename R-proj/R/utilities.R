modifyMat <- function(A){
  
  b=A[,1]
  A2=A[,-c(1)]
  retList=list("matrix"=A2, "vector"=b)
  return(retList)
  
}

#char2numeric <- function(r,){
  
#}

ineToMatrix <- function(P){
  
  #nrows=dim(P)[1]
  r=as.character(P[3,1])
  #print(r)
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
  #print(beg)
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
  #m=as.double(substr(r, start=2, stop=3))
  #d=as.double(substr(r, start=5, stop=6))
  A=rep(0,d)
  A[1]=m
  A[2]=d
  #print(m)
  #print(d)
  newrow=rep(0,d)
  #print(dim(P)[1]-2)
  #print(x[dim(P)[1]-2,1])
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
    #print(newrow)
    A=rbind(A,newrow)
    newrow=rep(0,d)
  }
  A=matrix(A,ncol=dim(A)[2])
  
  return(A)
}

testRvolEsti <- function(){
  path=getwd()
  path=paste0(substr(path, start=1, stop=nchar(path)-7),'/test/test_data/')
  print(path)
  listofexamples=list.files(path)
  
  for(i in 1:length(listofexamples)){
    x=read.csv(paste0(path,listofexamples[i]))
    print(listofexamples[i])
    A=ineToMatrix(x)
    tim=proc.time()
    VolEsti(list("matrix"=A))
    tim=proc.time()-tim
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  #return(A)
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
  verbose=FALSE
  if(!is.null(Inputs$verbose)){
    if(Inputs$verbose){
      verbose=TRUE
    }else{
      verbose=FALSE
    }
  }
  #print(A)
  #print(b)
  A=matrix(cbind(b,A),ncol=dim(A)[2]+1)
  #print(A)
  A=matrix(rbind(r,A),ncol=dim(A)[2])
  #return(list("matrix"=A,"vector"=b,"cheb"=xc))
  return(vol_R(A,10,1,xc,verbose))
  
}
