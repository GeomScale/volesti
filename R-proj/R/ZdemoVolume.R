

ZdemoVolume <- function(algo){
  
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
  print('--------1st test [2-dimensional]--------')
  print('----------------------------------------')
  cat('\n')
  print('Zonotope_2_4')
  ZonoMat = GenZonotope(2, 4)
  runtest(ZonoMat, 'Zonotope_2_4', tol, num_of_exps, algo)
  
  print('Zonotope_2_8')
  ZonoMat = GenZonotope(2, 8)
  runtest(ZonoMat, 'Zonotope_2_8', tol, num_of_exps, algo)
  
  print('----------------------------------------')
  print('-------2nd test [4-dimensional]---------')
  print('----------------------------------------')
  cat('\n')
  
  print('Zonotope_4_8')
  ZonoMat = GenZonotope(4, 8)
  runtest(ZonoMat, 'Zonotope_4_8', tol, num_of_exps, algo)
  
  print('Zonotope_4_10')
  ZonoMat = GenZonotope(4, 10)
  runtest(ZonoMat, 'Zonotope_4_10', tol, num_of_exps, algo)
  
  print('----------------------------------------')
  print('--------3rd test [5-dimensional]---------')
  print('----------------------------------------')
  cat('\n')
  print('Zonotope_5_10')
  ZonoMat = GenZonotope(5, 10)
  runtest(ZonoMat, 'Zonotope_5_10', tol, num_of_exps, algo)
  
}

runtest <- function(Mat, name_string, tol, num_of_exps, algo){
  
  exactvol = ExactZonoVol(Mat)
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + volume(list("matrix"=Mat, "Zonotope"=TRUE))
    } else {
      vol = vol + volume(list("matrix"=Mat, "Zonotope"=TRUE, "CG"=TRUE, "error"=0.2))
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

