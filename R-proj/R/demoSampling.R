#' Run some sampling experiments.
#'
#' Use uniform or spherical gaussian to sample from some convex H-polytopes, i.e. cubes, simplices, skinny cubes, cross polytopes and birkhoff polytopes.
#' We use the default values, i.e. \eqn{walk length = \lfloor 10+dimension/10\rfloor}, \eqn{N = 100}, Cordinate Directions HnR, \eqn{variance = 1}.
#' 
#' @param uniform The string "uniform" to choose uniform as the target distribution.
#' @param gaussian The string "gaussian" to choose spherical gaussian as the target distribution.
#' 
#' @return Print the computed volumes and the error. If the test fails a message is printed.
#' @examples
#' # choose uniform distribution
#' demoSampling("uniform")
#' # choose spherical gaussian distribution
#' demoSampling("gaussian")
demoSampling <- function(distribution){
  
  if(distribution!="uniform" & distribution!="gaussian"){
    print('choose between uniform and gaussian for the distribution.')
    return()
  }
  
  path = getwd()
  path = paste0(substr(path, start=1, stop=nchar(path) - 7), '/data/')
  
  print('----------------------------------------')
  print('------------1st test [cubes]------------')
  print('----------------------------------------')
  cat('\n')
  PolyList = GenCube(10, 'H')
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-cube10', distribution)
  
  PolyList = GenCube(20, 'H')
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-cube20', distribution)
  
  PolyMat = GenCube(5, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cube5', distribution)
  
  PolyMat = GenCube(10, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cube10', distribution)

  print('----------------------------------------')
  print('------2nd test [cross_polytopes]--------')
  print('----------------------------------------')
  cat('\n')
  PolyList = GenCross(10, 'H')
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-cross10', distribution)
  
  PolyMat = GenCross(10, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cross10', distribution)
  
  print('----------------------------------------')
  print('----------3rd test [birkhoff]-----------')
  print('----------------------------------------')
  cat('\n')
  A = ineToMatrix(read.csv(paste0(path,'birk3.ine')))
  x = modifyMat(A)
  runsample(x$matrix, x$vector, FALSE, FALSE, 'H-birk3', distribution)
  
  A = ineToMatrix(read.csv(paste0(path,'birk4.ine')))
  x = modifyMat(A)
  runsample(x$matrix, x$vector, FALSE, FALSE, 'H-birk4', distribution)
  
  A = ineToMatrix(read.csv(paste0(path,'birk5.ine')))
  x = modifyMat(A)
  runsample(x$matrix, x$vector, FALSE, FALSE, 'H-birk5', distribution)
  
  A = ineToMatrix(read.csv(paste0(path,'birk6.ine')))
  x = modifyMat(A)
  runsample(x$matrix, x$vector, FALSE, FALSE, 'H-birk6', distribution)
  
  print('----------------------------------------')
  print('--------4th test [prod_simplex]---------')
  print('----------------------------------------')
  cat('\n')
  PolyList = GenProdSimplex(5)
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-prod_simplex_5_5', distribution)
  
  PolyList = GenProdSimplex(10)
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-prod_simplex_10_10', distribution)
  
  PolyList = GenProdSimplex(15)
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-prod_simplex_15_15', distribution)
  
  PolyList = GenProdSimplex(20)
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-prod_simplex_20_20', distribution)
  
  print('----------------------------------------')
  print('-----------5th test [simplex]-----------')
  print('----------------------------------------')
  cat('\n')
  PolyList = GenSimplex(10, 'H')
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-prod_simplex10', distribution)
  
  PolyList = GenSimplex(20, 'H')
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-simplex20', distribution)
  
  PolyMat = GenSimplex(10, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-simplex10', distribution)
  
  PolyMat = GenSimplex(20, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-simplex20', distribution)
  
  print('----------------------------------------')
  print('--------6th test [skinny_cubes]---------')
  print('----------------------------------------')
  cat('\n')
  PolyList = GenSkinnyCube(10)
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-skinny_cube10', distribution)
  
  PolyList = GenSkinnyCube(20)
  runsample(PolyList$matrix, PolyList$vector, FALSE, FALSE, 'H-skinny_cube20', distribution)
  
  print('----------------------------------------')
  print('----------7th test [Zonotopes]----------')
  print('----------------------------------------')
  cat('\n')
  zonotope = GenZonotope(4, 8)
  runsample(zonotope, c(0,0), FALSE, TRUE, 'zonotope_4_8', distribution)
  
  zonotope = GenZonotope(5, 10)
  runsample(zonotope, c(0,0), FALSE, TRUE, 'zonotope_5_10', distribution)
  
}


runsample <- function(Mat, vector, Vpoly, Zono, name_string, dist){
  if (dist == "uniform") {
    if (Zono) {
      p = sample_points(list("matrix"=Mat, "Zonotope"=TRUE))
    } else if (Vpoly) {
      p = sample_points(list("matrix"=Mat, "Vpoly"=TRUE))
    } else {
      p = sample_points(list("matrix"=Mat, "vector"=vector))
    }
  } else {
    if (Zono) {
      p = sample_points(list("matrix"=Mat, "Zonotope"=TRUE, "gaussian"=TRUE))
    } else if (Vpoly) {
      p = sample_points(list("matrix"=Mat, "Vpoly"=TRUE, "gaussian"=TRUE))
    } else {
      p = sample_points(list("matrix"=Mat, "vector"=vector, "gaussian"=TRUE))
    }
  }
  if (length(p[is.nan(p)])>0 | length(p[is.infinite(p)])>0) {
    print(paste0('Test FAILED!! [', name_string, ']  There are NaN or Inf values in the coordinates of your points!'))
  } else {
    print(paste0('Test PASSED!! [', name_string, ']'))
  }
  cat('\n')
}
