#' function to get a ine file and returns a numerical matrix A.
#'
#' This function takes an ine file as a string (using read.csv()) and returns a numerical matrix A in ine format for function volume (see \eqn{volume} function examples).
#'
#' @param P It is in format, read.cs('path/to/file.ine'). The ine file desrcibes a H-polytope.
#' @return The numerical matrix in ine format.
#' @examples
#' #give the path to cube40.ine
#' A = ineToMatrix(read.csv('path/to/data/cube40.ine'))
ineToMatrix <- function(P){
  
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
  A = matrix(A, ncol=dim(A)[2])
  
  return(A)
}
