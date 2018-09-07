#' function to get a ine file and returns a numerical matrix A.
#'
#' This function takes the path for an ine or an ext file and returns the corresponding numerical matrix and vector that are compatible with volesti package's functions.
#'
#' @param path A string that containes the path to an ine or a ext file. The ine file desrcibes a H-polytope and ext file describes a V-polytope or a zonotope.
#' @return If the path corresponds to an ine file then the return value is a list that containes elements "A" and "b", i.e. the numerical \eqn{m\times d} matrix \eqn{A} and the numerical \eqn{m}-dimensional vector \eqn{b}, defining H-polytope \eqn{P}, s.t.:  \eqn{Ax\leq b}. If it corresponds to an ext file (V-polytopes or zonotopes) then the return value is a \eqn{m\times d} matrix that containes row-wise the vertices or the segments respectively.
#' @examples
#' # give the path to birk4.ine
#' ListPoly = fileToMatrix(path/to/data/birk4.ine)
fileToMatrix <- function(path){
  
  ineorext=substr(path, start = nchar(path) - 2, stop = nchar(path))
  P = read.csv(path)
  r = as.character(P[3,1])
  count_sp = 1
  str = ""
  beg = 0
  for (j in 1:nchar(r)) {
    if (substr(r, start=j, stop=j) == " ") {
      beg = beg + 1
    } else {
      break
    }
  }
  for (i in seq(from= beg + 1, to=nchar(r), by=1)) {
    if (substr(r, start=i, stop=i) == " ") {
      if (count_sp == 1) {
        m = as.numeric(str)
        str = ""
        count_sp = count_sp + 1
      } else {
        d = as.numeric(str)
        str = ""
        break
      }
    } else {
      str = paste0(str, substr(r, start=i, stop=i))
    }
  }
  A = rep(0,d)
  A[1] = m
  A[2] = d
  newrow = rep(0,d)
  for (i in 4:(dim(P)[1] - 2)) {
    r = P[i,1]
    r = as.character(r)
    str = ""
    count = 1
    beg = 0
    for (j in 1:nchar(r)) {
      if(substr(r, start=j, stop=j)==" "){
        beg = beg + 1
      } else {
        break
      }
    }
    sp_bef = FALSE
    for (j in seq(from=beg + 1, to=nchar(r), by=1)) {
      if (substr(r, start=j, stop=j) == " "){
        if (sp_bef) {
          next
        }
        sp_bef = TRUE
        newrow[count] = as.numeric(str)
        str = ""
        count = count + 1
      } else {
        str = paste0(str, substr(r, start=j, stop=j))
        sp_bef = FALSE
        if (j == nchar(r)) {
          newrow[count] = as.numeric(str)
        }
      }
    }
    A = rbind(A,newrow)
    newrow = rep(0,d)
  }
  A = matrix(A, ncol=dim(A)[2]) # now matrix A is in ine or ext format
  
  # remove first row
  A = A[-c(1),]
  
  # first column is the vector b
  b = A[,1]
  
  # remove first column
  A2 = A[,-c(1)]
  
  if(ineorext=="ine") {
    # return final list if it is ine file
    retList = list("A"=A2, "b"=b)
    return(retList)
  } else if(ineorext == "ext") {
    return(A2)
  } else {
    print('Only ine or et files can be handled by this function!')
    return(-1)
  }
}
