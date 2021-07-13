#' Read a SDPA format file
#'
#' Read a SDPA format file and return a spectrahedron (an object of class Spectrahedron) which is defined by
#' the linear matrix inequality in the input file, and the objective function.
#'
#' @param path Name of the input file
#'
#' @return A list with two named items: an item "matrices" which is an object of class Spectrahedron and an vector "objFunction"
#'
#' @examples
#' path = system.file('extdata', package = 'volesti')
#' l = read_sdpa_format_file(paste0(path,'/sdpa_n2m3.txt'))
#' Spectrahedron = l$spectrahedron
#' objFunction = l$objFunction
#' @export
#' @useDynLib volesti, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp loadModule
#' @importFrom "methods" "new"
#' @exportPattern "^[[:alpha:]]+"
read_sdpa_format_file <- function(path){
    l = load_sdpa_format_file(path)
    S = Spectrahedron(matrices = l$matrices)

    return(list("spectrahedron"=S, "objFunction"= l$objFunction))
}
