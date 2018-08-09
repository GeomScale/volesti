#' Run some sampling experiments.
#'
#' Use uniform or spherical gaussian to sample from some convex H-polytopes, i.e. cubes, simplices, skinny_cubes, cross polytopes, birkhoff polytopes.
#' We use the default values, i.e. \eqn{walk length = \lfloor 10+dimension/10\rfloor}, \eqn{N = 100}, Cordinate Directions HnR, \eqn{variance = 1}.
#' 
#' @param uniform The string "uniform" to choose uniform as the target distribution.
#' @param gaussian The string "gaussian" to choose spherical gaussian as the target distribution.
#' 
#' @return Print the computed volumes and the error. If the test fails a message is printed.
#' @examples
#' #choose uniform distribution
#' demoSampling("uniform")
#' #choose spherical gaussian distribution
#' demoSampling("gaussian")
demoSampling <- function(distribution){
  if(distribution!="uniform" & distribution!="gaussian"){
    print('choose between uniform and gaussian for the distribution.')
    return()
  }
  
  path=getwd()
  path=paste0(substr(path, start=1, stop=nchar(path)-7),'/data/')
  
  print('----------------------------------------')
  print('------------1st test [cubes]------------')
  print('----------------------------------------')
  listofpoly=c('cube10.ine','cube20.ine','cube30.ine')
  runsample(path, listofpoly, distribution)

  print('----------------------------------------')
  print('------2nd test [cross_polytopes]--------')
  print('----------------------------------------')
  listofpoly=c('cross_10.ine')
  runsample(path, listofpoly, distribution)
  
  print('----------------------------------------')
  print('----------3rd test [birkhoff]-----------')
  print('----------------------------------------')
  listofpoly=c('birk3.ine','birk4.ine','birk5.ine','birk6.ine')
  runsample(path, listofpoly, distribution)
  
  print('----------------------------------------')
  print('--------4th test [prod_simplex]---------')
  print('----------------------------------------')
  cat('\n')
  listofpoly=c('prod_simplex_5_5.ine','prod_simplex_10_10.ine','prod_simplex_15_15.ine','prod_simplex_20_20.ine')
  runsample(path, listofpoly, distribution)
  
  print('----------------------------------------')
  print('-----------5th test [simplex]-----------')
  print('----------------------------------------')
  listofpoly=c('simplex10.ine','simplex20.ine','simplex30.ine','simplex40.ine','simplex50.ine')
  runsample(path, listofpoly, distribution)
  
  print('----------------------------------------')
  print('--------6th test [skinny_cubes]---------')
  print('----------------------------------------')
  listofpoly=c('skinny_cube10.ine','skinny_cube20.ine')
  runsample(path, listofpoly, distribution)
}


runsample <- function(path,listofpoly,dist){
  for(i in 1:length(listofpoly)){
    if (dist == "uniform"){
      p = sample_points(list("path"=paste0(path,listofpoly[i])))
    } else {
      p = sample_points(list("path"=paste0(path,listofpoly[i]), "gaussian"=TRUE))
    }
    if (length(p[is.nan(p)])>0 | length(p[is.infinite(p)])>0){
      print(paste0('Test FAILED!! [',listofpoly[i],']  There are NaN or Inf values in the coordinates of your points!'))
    } else {
      print(paste0('Test PASSED!! [',listofpoly[i],']'))
    }
  }
  cat('\n')
}