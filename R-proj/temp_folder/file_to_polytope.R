#' function to get an ine or an ext file and returns the corresponding polytope
#'
#' For an ".ine" file it generates the corresponding H-polytope. For an ".ext" file it generates the corresponding V-polytope or zonotope.
#' For more details on those file formats see \url{https://github.com/GeomScale/volume_approximation/blob/develop/doc/cpp_interface.md#polytope-input}.
#'
#' @param path A string that containes the path to an ine or a ext file. The ine file desrcibes a H-polytope and ext file describes a V-polytope or a zonotope.
#' @param zonotope A boolean parameter. It has to be TRUE when the path leads to an .ext file that describes a zonotope.
#'
#' @return A polytope class. If the path corresponds to an ine file then the return value represents a H-polytope. If it corresponds to an ext file the return value represents a V-polytope (default choice) or a zonotope if the second argument is TRUE.
#'
#' @export
#' @useDynLib volesti, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp loadModule
#' @importFrom "utils" "read.csv"
#' @importFrom "stats" "cov"
#' @importFrom "methods" "new"
#' @exportPattern "^[[:alpha:]]+"
file_to_polytope <- function(path, zonotope = FALSE){
  
  ineorext=substr(path, start = nchar(path) - 2, stop = nchar(path))
  if(ineorext!="ine" && ineorext!="ext") {
    stop("Only ine or ext files can be handled by this function!")
  }
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
    P = Hpolytope$new(-A2,b)
  } else {
    if(!missing(zonotope)){
      if(zonotope) {
        P = Zonotope$new(A2)
        return(P)
      }
    }
    P = Vpolytope$new(A2)
  }
  return(P)
}
