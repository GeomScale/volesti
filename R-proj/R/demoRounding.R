#' Run rounding and rotating tests.
#' 
#' Choose volume algorithm between CoolingGaussian and SequenceOfBalls and run rounding tests for some skinny cubes. In the first test we apply a random rotation as well before the rounding. For we run 10 experiments for SequenceOfBalls and 20 for CoolingGaussian.
#' 
#' @param CG The string "CG" to choose CoolingGaussian algorithm
#' @param SOB The string "SOB" to choose SequenceOfBalls algorithm
#' @return Print the computed volume and print a failure message if the error is larger than the expected.
#' 
#' @examples 
#' #run tests for volesti algorithm
#' demoRounding("volesti")
#' 
#' #run tests for CV algorithm
#' demoRounding("CV")
demoRounding <- function(algo){
  
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
  
  print('--------------------------------------------')
  print('------1st test [rotated_skinny_cubes]-------')
  print('--------------------------------------------')
  cat('\n')
  PolyList = GenSkinnyCube(10)
  testRound(PolyList$matrix, PolyList$vector, 102400, 0.1, 'H-skinny_cube10', num_of_exps, algo, TRUE)
  
  PolyList = GenSkinnyCube(20)
  testRound(PolyList$matrix, PolyList$vector, 104857600, 0.3, 'H-skinny_cube20', num_of_exps, algo, TRUE)
  
  print('--------------------------------------------')
  print('------2nd test [skinny_cubes]-------')
  print('--------------------------------------------')
  cat('\n')
  PolyList = GenSkinnyCube(10)
  testRound(PolyList$matrix, PolyList$vector, 102400, 0.1, 'H-skinny_cube10', num_of_exps, algo, FALSE)
  
  PolyList = GenSkinnyCube(20)
  testRound(PolyList$matrix, PolyList$vector, 104857600, 0.3, 'H-skinny_cube20', num_of_exps, algo, FALSE)
  
}

testRound <- function(Mat, vector, exactvol, tol, name_string, num_of_exps, algo, rotation){
  
  #for (i in 1:length(listofpoly)) {
    if (rotation) {
      listHpoly = rand_rotate(list("matrix"=Mat, "vector"=vector))
      listHpoly = round_polytope(list("matrix"=listHpoly$matrix, "vector"=listHpoly$vector))
    } else {
      listHpoly = round_polytope(list("matrix"=Mat, "vector"=vector))
    }
    vol = 0
    for (j in 1:num_of_exps) {
      if (algo == "SOB") {
        vol = vol + listHpoly$round_value * volume(list("matrix"=listHpoly$matrix, "vector"=listHpoly$vector))
      } else {
        vol = vol + listHpoly$round_value * volume(list("matrix"=listHpoly$matrix, "vector"=listHpoly$vector, "CG"=TRUE, "error"=0.1))
      }
    }
    vol = vol / num_of_exps
    print(paste0('volume approximation of ', name_string,' = ',vol))
    print(paste0('exact volume of ', name_string,' = ',exactvol))
    error = abs(vol - exactvol) / exactvol
    print(paste0('error = ',error))
    if (error >= tol){
      print(paste0('TEST FAILED!! ', error, ' > ', tol))
    }
    cat('\n')
  #}
}
