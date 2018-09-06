#' Run some volume approximation experiments for H-polytopes.
#'
#' Choose between SequenceOfBalls and CoolingGaussian algorithm to approximate the volume of some cubes, simplices, skinny_cubes, cross polytopes and birkhoff polytopes in H-representation.
#' For each polytope we run \eqn{10} volume experiments for SequenceOfBalls and \eqn{20} for CoolingGaussian and we consider the mean value as the volume approximation. We demand \eqn{error = 0.1} for the most of them. For all the other parameters we use the default values for both algorithms.
#' 
#' @param CG The string "CG" to choose CoolingGaussian algorithm.
#' @param SOB The string "SOB" to choose SequenceOfBalls algorithm.
#' 
#' @return Print the computed volumes and the error. If the test fails a message is printed.
#' @examples
#' # test SequenceOfBalls
#' HdemoVolume("SOB")
#' # test CoolingGausian
#' HdemoVolume("CG")
HdemoVolume <- function(algo){
  
  if(algo!="CG" & algo!="SOB"){
    print('choose between CV and volesti.')
    return()
  }
  if (algo == "CG") {
    num_of_exps = 20
  } else {
    num_of_exps = 10
  }
  path = getwd()
  path = paste0(substr(path, start=1, stop=nchar(path) - 7), '/data/')
  
  print('----------------------------------------')
  print('------------1st test [cubes]------------')
  print('----------------------------------------')
  cat('\n')
  print('H-cube10')
  PolyList = GenCube(10, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-cube10', 1024, 0.1, num_of_exps, algo)
  
  print('H-cube20')
  PolyList = GenCube(20, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-cube20', 1048576, 0.1, num_of_exps, algo)
  
  print('H-cube30')
  PolyList = GenCube(30, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-cube30', 1073742000, 0.2, num_of_exps, algo)
  
  print('----------------------------------------')
  print('------2nd test [cross_polytopes]--------')
  print('----------------------------------------')
  cat('\n')
  
  print('H-cross10')
  PolyList = GenCross(10, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-cross10', 0.0002821869, 0.1, num_of_exps, algo)
  
  print('----------------------------------------')
  print('----------3rd test [birkhoff]-----------')
  print('----------------------------------------')
  cat('\n')
  print('H-birk3')
  A = ineToMatrix(read.csv(paste0(path,'birk3.ine')))
  x = modifyMat(A)
  runtest(x$matrix, x$vector, 'H-birk3', 0.125, 0.1, num_of_exps, algo)
  
  print('H-birk4')
  A = ineToMatrix(read.csv(paste0(path,'birk4.ine')))
  x = modifyMat(A)
  runtest(x$matrix, x$vector, 'H-birk4', 0.000970018, 0.2, num_of_exps, algo)
  
  print('H-birk5')
  A = ineToMatrix(read.csv(paste0(path,'birk5.ine')))
  x = modifyMat(A)
  runtest(x$matrix, x$vector, 'H-birk5', 0.000000225, 0.2, num_of_exps, algo)
  
  print('H-birk6')
  A = ineToMatrix(read.csv(paste0(path,'birk6.ine')))
  x = modifyMat(A)
  
  runtest(x$matrix, x$vector, 'H-birk6', 0.0000000000009455459196, 0.5, num_of_exps, algo)
  
  print('----------------------------------------')
  print('--------4th test [prod_simplex]---------')
  print('----------------------------------------')
  cat('\n')
  print('H-prod_simplex_5_5')
  PolyList = GenProdSimplex(5)
  runtest(PolyList$matrix, PolyList$vector, 'H-prod_simplex_5_5', (1/prod(1:5))^2, 0.1, num_of_exps, algo)
  
  print('H-prod_simplex_10_10')
  PolyList = GenProdSimplex(10)
  runtest(PolyList$matrix, PolyList$vector, 'H-prod_simplex_10_10', (1/prod(1:10))^2, 0.1, num_of_exps, algo)
  
  print('H-prod_simplex_15_15')
  PolyList = GenProdSimplex(15)
  runtest(PolyList$matrix, PolyList$vector, 'H-prod_simplex_15_15', (1/prod(1:15))^2, 0.1, num_of_exps, algo)
  
  print('H-prod_simplex_20_20')
  PolyList = GenProdSimplex(20)
  runtest(PolyList$matrix, PolyList$vector, 'H-prod_simplex_20_20', (1/prod(1:20))^2, 0.1, num_of_exps, algo)
  
  print('----------------------------------------')
  print('--------5th test [simplex]---------')
  print('----------------------------------------')
  cat('\n')
  print('H-simplex10')
  PolyList = GenSimplex(10, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-simplex10', (1/prod(1:10)), 0.1, num_of_exps, algo)
  
  print('H-simplex20')
  PolyList = GenSimplex(20, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-simplex20', (1/prod(1:20)), 0.1, num_of_exps, algo)
  
  print('H-simplex30')
  PolyList = GenSimplex(30, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-simplex30', (1/prod(1:30)), 0.1, num_of_exps, algo)
  
  print('H-simplex40')
  PolyList = GenSimplex(40, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-simplex40', (1/prod(1:40)), 0.1, num_of_exps, algo)
  
  print('H-simplex50')
  PolyList = GenSimplex(50, 'H')
  runtest(PolyList$matrix, PolyList$vector, 'H-simplex50', (1/prod(1:50)), 0.1, num_of_exps, algo)
  
  if(algo=="SOB"){
    print('----------------------------------------')
    print('--------6th test [skinny_cubes]---------')
    print('----------------------------------------')
    cat('\n')
    print('H-skinny_cube10')
    PolyList = GenSkinnyCube(10)
    runtest(PolyList$matrix, PolyList$vector, 'H-skinny_cube10', 102400, 0.1, num_of_exps, algo)
    
    print('H-skinny_cube20')
    PolyList = GenSkinnyCube(20)
    runtest(PolyList$matrix, PolyList$vector, 'H-skinny_cube20', 104857600, 0.1, num_of_exps, algo)
  }
  
}

runtest <- function(Mat, vector, name_string, exactvol, tol, num_of_exps, algo){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + volume(list("matrix"=Mat, "vector"=vector))
    } else {
      vol = vol + volume(list("matrix"=Mat, "vector"=vector, "CG"=TRUE, "error"=0.1))
    }
  }
  vol = vol / num_of_exps
  print(paste0('volume approximation of ',name_string,' = ',vol))
  print(paste0('exact volume of ',name_string,' = ',exactvol))
  error = abs(vol - exactvol) / exactvol
  print(paste0('error = ',error))
  if (error >= tol){
    print(paste0('TEST FAILED!! ', error, ' > ', tol))
  }
  cat('\n')
}
