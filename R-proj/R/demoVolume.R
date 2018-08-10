#' Run some volume approxiamtion experiments.
#'
#' Run volesti or CV algorithm to approximate the volume of some cubes, simplices, skinny_cubes, cross polytopes, birkhoff polytopes.
#' We run \eqn{10} experiments for volesti and \eqn{20} for CV, for each polytope and we consider the mean as the computed volume. We demand \eqn{error = 0.1}. For all the other parameters use the default values for both algorithms.
#' 
#' @param CV The string "CV" to choose CV algorithm.
#' @param volesti The string "volesti" to choose volesti algorithm.
#' 
#' @return Print the computed volumes and the error. If the test fails a message is printed.
#' @examples
#' #test volesti
#' demoVolume("volesti")
#' #test CV
#' demoVolume("CV")
demoVolume <- function(algo){
  
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
  
  print('----------------------------------------')
  print('------------1st test [cubes]------------')
  print('----------------------------------------')
  cat('\n')
  listofpoly = c('cube10.ine', 'cube20.ine', 'cube30.ine')
  exactvols=c(1024, 1048576, 1073742000)
  tols = c(0.1, 0.1, 0.2)
  runtest(path, listofpoly, exactvols, tols, num_of_exps, algo)
  
  
  print('----------------------------------------')
  print('------2nd test [cross_polytopes]--------')
  print('----------------------------------------')
  cat('\n')
  listofpoly=c('cross_10.ine')
  exactvols=c(0.0002821869)
  tols=c(0.1)
  runtest(path, listofpoly, exactvols, tols, num_of_exps, algo)
  
  
  print('----------------------------------------')
  print('----------3rd test [birkhoff]-----------')
  print('----------------------------------------')
  cat('\n')
  listofpoly=c('birk3.ine', 'birk4.ine', 'birk5.ine', 'birk6.ine')
  exactvols=c(0.125, 0.000970018, 0.000000225, 0.0000000000009455459196)
  tols=c(0.1, 0.2, 0.2, 0.5)
  runtest(path, listofpoly, exactvols, tols, num_of_exps, algo)
  
  
  print('----------------------------------------')
  print('--------4th test [prod_simplex]---------')
  print('----------------------------------------')
  cat('\n')
  listofpoly=c('prod_simplex_5_5.ine', 'prod_simplex_10_10.ine', 'prod_simplex_15_15.ine', 'prod_simplex_20_20.ine')
  exactvols=c( (1/prod(1:5))^2, (1/prod(1:10))^2, (1/prod(1:15))^2, (1/prod(1:20))^2)
  tols=c(0.1, 0.1, 0.1, 0.1)
  runtest(path, listofpoly, exactvols, tols, num_of_exps, algo)
  
  
  print('----------------------------------------')
  print('--------5th test [simplex]---------')
  print('----------------------------------------')
  cat('\n')
  listofpoly=c('simplex10.ine', 'simplex20.ine', 'simplex30.ine', 'simplex40.ine', 'simplex50.ine')
  exactvols=c( 1/prod(1:10), 1/prod(1:20), 1/prod(1:30), 1/prod(1:40), 1/prod(1:50))
  tols=c(0.1, 0.1, 0.1, 0.1, 0.1)
  runtest(path, listofpoly, exactvols, tols, num_of_exps, algo)
  
  
  if(algo=="volesti"){
    print('----------------------------------------')
    print('--------6th test [skinny_cubes]---------')
    print('----------------------------------------')
    cat('\n')
    listofpoly=c('skinny_cube10.ine', 'skinny_cube20.ine')
    exactvols=c(102400, 104857600)
    tols=c(0.1, 0.1)
    runtest(path, listofpoly, exactvols, tols, num_of_exps, algo)
  }
  
}



runtest <- function(path, listofpoly, exactvols, tols, num_of_exps, algo){
  
  for (i in 1:length(listofpoly)) {
    vol = 0
    for (j in 1:num_of_exps) {
      if (algo == "volesti") {
        vol = vol + volume(list("path"=paste0(path,listofpoly[i])))
      } else {
        vol = vol + volume(list("path"=paste0(path,listofpoly[i]), "CV"=TRUE, "error"=0.1))
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
