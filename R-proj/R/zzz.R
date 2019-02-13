
## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}

loadModule("yada", TRUE)

## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
#' Classes to construct convex polytopes (H, V or zonotopes)
#'
#' C++ classes exposed in R by rcpp_module
#' @usage Hpolytope$new(A,b) or Vpolytope$new(V) or Zonotope$new(G)
#' 
#' @examples
#' # Create a 2-d unit simplex in H-representation
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = Hpolytope$new(A,b)
#'
#' # Create a 3-d cube in V-representation
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' P = Vpolytope$new(V)
#'
#' # Create a 2-d zonotope with 4 generators
#' G = matrix(c(1,0,0,1,-0.73,0.67,-0.25,0.96), ncol = 2, nrow = 4, byrow = TRUE)
#' P = Zonotope$new(G)
Polytopes <- function(A,b){
    #empty function. Only for Rd file for polytope modules
}


