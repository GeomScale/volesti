#' Run some volume approximation experiments for V-polytopes.
#'
#' Choose between SequenceOfBalls and CoolingGaussian algorithm to approximate the volume of some cubes, simplices and cross polytopes in V-representation.
#' For each polytope we run \eqn{10} volume experiments and we consider the mean value as the volume approximation. For SOB algorithm we demand \eqn{error = 0.1} and for CG algorithm we demand \eqn{error = 0.2}.
#' 
#' @param CG The string "CG" to choose CoolingGaussian algorithm.
#' @param SOB The string "SOB" to choose SequenceOfBalls algorithm.
#' 
#' @return Print the computed volumes and the error. If the test fails a message is printed.
#' @examples
#' # test SequenceOfBalls
#' VdemoVolume("SOB")
#' # test CoolingGausian
#' VdemoVolume("CG")
VdemoVolume <- function(algo){
  
  if(algo!="CG" & algo!="SOB"){
    print('choose between CV and volesti.')
    return()
  }
  
  num_of_exps = 10
  if (algo == 'CG') {
    tol = 0.2
  } else {
    tol = 0.1
  }
  
  print('----------------------------------------')
  print('------------1st test [cubes]------------')
  print('----------------------------------------')
  cat('\n')
  print('V-cube3')
  PolyMat = GenCube(3, 'V')
  Vruntest(PolyMat, 'V-cube3', 8, tol, num_of_exps, algo)
  
  print('V-cube4')
  PolyMat = GenCube(4, 'V')
  Vruntest(PolyMat, 'V-cube4', 16, tol, num_of_exps, algo)
  
  print('----------------------------------------')
  print('------2nd test [cross_polytopes]--------')
  print('----------------------------------------')
  cat('\n')
  
  print('V-cross5')
  PolyMat = GenCross(5, 'V')
  Vruntest(PolyMat, 'V-cross5', 0.2666667, tol, num_of_exps, algo)
  
  print('V-cross7')
  PolyMat = GenCross(7, 'V')
  Vruntest(PolyMat, 'V-cross7', 0.02539683, tol, num_of_exps, algo)
  
  print('----------------------------------------')
  print('--------3rd test [simplex]---------')
  print('----------------------------------------')
  cat('\n')
  print('V-simplex5')
  PolyMat = GenSimplex(5, 'V')
  Vruntest(PolyMat, 'V-simplex5', 1/prod(1:5), tol, num_of_exps, algo)
  
  print('V-simplex7')
  PolyMat = GenSimplex(7, 'V')
  Vruntest(PolyMat, 'V-simplex7', 1/prod(1:7), tol, num_of_exps, algo)
  

  Vruntest <- function(Mat, name_string, exactvol, tol, num_of_exps, algo){
  
    vol = 0
    for (j in 1:num_of_exps) {
      if (algo == "SOB") {
        vol = vol + volume(list("matrix"=Mat, "Vpoly"=TRUE))
      } else {
        vol = vol + volume(list("matrix"=Mat, "Vpoly"=TRUE, "CG"=TRUE, "error"=0.2))
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

}
