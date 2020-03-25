#' Read a SDPA format file
#'
#' Read a SDPA format file and return a spectrahedron (an object of class Spectrahedron) which is defined by
#' the linear matrix inequality in the input file, and the objective function.
#'
#' @field path Name of the input file
#'
#' @return A list with two named items: an item "matrices" which is an object of class Spectrahedron and an vector "objFunction"
#'
#' @examples
#' l = loadSdpaFormatFile("input.txt")
#' @export
#' @useDynLib volesti, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp loadModule
#' @importFrom "methods" "new"
#' @exportPattern "^[[:alpha:]]+"
readSdpaFormatFile <- function(path){
    l = loadSdpaFormatFile(path)
    S = Spectrahedron$new(l$matrices)

    return(list("spectrahedron"=S, "objFunction"= l$objFunction))
}
