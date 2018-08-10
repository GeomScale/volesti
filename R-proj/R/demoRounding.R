#' Run rounding and rotating tests.
#' 
#' Choose volume algorithm between CV and volesti and run rounding tests for some skinny cubes. In the first test we apply a random rotation as well before the rounding. For we run 10 experiments for volesti and 20 for CV.
#' 
#' @param CV The string "CV" to choose CV algorithm
#' @param volesti The string "volesti" to choose volesti algorithm
#' @return Print the computed volume and print a failure message if the error is larger than the expected.
#' 
#' @examples 
#' #run tests for volesti algorithm
#' demoRounding("volesti")
#' 
#' #run tests for CV algorithm
#' demoRounding("CV")
demoRounding <- function(algo){
  
  if(algo!="CV" & algo!="volesti"){
    print('choose between CV and volesti.')
    return()
  }
  if (algo == "CV") {
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
  listofpoly=c('skinny_cube10.ine', 'skinny_cube20.ine')
  exactvols=c(102400, 104857600)
  tols=c(0.1, 0.3)
  
  testRound(path, listofpoly, exactvols, tols, num_of_exps, algo, TRUE)
  
  print('--------------------------------------------')
  print('------2nd test [skinny_cubes]-------')
  print('--------------------------------------------')
  cat('\n')
  listofpoly=c('skinny_cube10.ine', 'skinny_cube20.ine')
  exactvols=c(102400, 104857600)
  tols=c(0.1, 0.1)
  testRound(path, listofpoly, exactvols, tols, num_of_exps, algo, FALSE)
  
}

testRound <- function(path, listofpoly, exactvols, tols, num_of_exps, algo, rotation){
  
  for (i in 1:length(listofpoly)) {
    if (rotation) {
      listHpoly = rand_rotate(list("path"=paste0(path,listofpoly[i])))
      listHpoly = round_polytope(list("matrix"=listHpoly$matrix, "vector"=listHpoly$vector))
    } else {
      listHpoly = round_polytope(list("path"=paste0(path,listofpoly[i])))
    }
    vol = 0
    for (j in 1:num_of_exps) {
      if (algo == "volesti") {
        vol = vol + listHpoly$round_value * volume(list("matrix"=listHpoly$matrix, "vector"=listHpoly$vector))
      } else {
        vol = vol + listHpoly$round_value * volume(list("matrix"=listHpoly$matrix, "vector"=listHpoly$vector, "CV"=TRUE, "error"=0.1))
      }
    }
    vol = vol / num_of_exps
    print(paste0('volume approximation of ',listofpoly[i],' = ',vol))
    print(paste0('exact volume of ',listofpoly[i],' = ',exactvols[i]))
    error = abs(vol - exactvols[i]) / exactvols[i]
    print(paste0('error = ',error))
    if (error >= tols[i]){
      print(paste0('TEST FAILED!! ', error, ' > ', tols[i]))
    }
    cat('\n')
  }
}
