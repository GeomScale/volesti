Hruntest <- function(Mat, vector, name_string, exactvol, tol, num_of_exps, alg){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (alg == "SOB") {
      vol = vol + volume(A=Mat, b=vector)
    } else {
      vol = vol + volume(A=Mat, b=vector, CG=TRUE, error=0.1)
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

library(volesti)

path = system.file('extdata', package = 'volesti')

for (i in 1:2) {
  
  if (i==1) {
    algo = 'CG'
    num_of_exps = 20
  } else {
    algo = 'SOB'
    num_of_exps = 10
  }

  print('----------------------------------------')
  print('------------1st test [cubes]------------')
  print('----------------------------------------')
  cat('\n')
  print('H-cube10')
  PolyList = GenCube(10, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-cube10', 1024, 0.1, num_of_exps, algo)

  print('H-cube20')
  PolyList = GenCube(20, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-cube20', 1048576, 0.1, num_of_exps, algo)

  print('H-cube30')
  PolyList = GenCube(30, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-cube30', 1073742000, 0.2, num_of_exps, algo)

  print('----------------------------------------')
  print('------2nd test [cross_polytopes]--------')
  print('----------------------------------------')
  cat('\n')

  print('H-cross10')
  PolyList = GenCross(10, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-cross10', 0.0002821869, 0.1, num_of_exps, algo)

  print('----------------------------------------')
  print('----------3rd test [birkhoff]-----------')
  print('----------------------------------------')
  cat('\n')
  print('H-birk3')
  ListPoly = fileToMatrix(paste0(path,'/birk3.ine'))
  Hruntest(ListPoly$A, ListPoly$b, 'H-birk3', 0.125, 0.1, num_of_exps, algo)

  print('H-birk4')
  ListPoly = fileToMatrix(paste0(path,'/birk4.ine'))
  Hruntest(ListPoly$A, ListPoly$b, 'H-birk4', 0.000970018, 0.2, num_of_exps, algo)

  print('H-birk5')
  ListPoly = fileToMatrix(paste0(path,'/birk5.ine'))
  Hruntest(ListPoly$A, ListPoly$b, 'H-birk5', 0.000000225, 0.2, num_of_exps, algo)

  print('H-birk6')
  ListPoly = fileToMatrix(paste0(path,'/birk6.ine'))
  Hruntest(ListPoly$A, ListPoly$b, 'H-birk6', 0.0000000000009455459196, 0.5, num_of_exps, algo)

  print('----------------------------------------')
  print('--------4th test [prod_simplex]---------')
  print('----------------------------------------')
  cat('\n')
  print('H-prod_simplex_5_5')
  PolyList = GenProdSimplex(5)
  Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_5_5', (1/prod(1:5))^2, 0.1, num_of_exps, algo)

  print('H-prod_simplex_10_10')
  PolyList = GenProdSimplex(10)
  Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_10_10', (1/prod(1:10))^2, 0.1, num_of_exps, algo)

  print('H-prod_simplex_15_15')
  PolyList = GenProdSimplex(15)
  Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_15_15', (1/prod(1:15))^2, 0.1, num_of_exps, algo)

  print('H-prod_simplex_20_20')
  PolyList = GenProdSimplex(20)
  Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_20_20', (1/prod(1:20))^2, 0.1, num_of_exps, algo)

  print('----------------------------------------')
  print('--------5th test [simplex]---------')
  print('----------------------------------------')
  cat('\n')
  print('H-simplex10')
  PolyList = GenSimplex(10, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-simplex10', (1/prod(1:10)), 0.1, num_of_exps, algo)

  print('H-simplex20')
  PolyList = GenSimplex(20, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-simplex20', (1/prod(1:20)), 0.1, num_of_exps, algo)

  print('H-simplex30')
  PolyList = GenSimplex(30, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-simplex30', (1/prod(1:30)), 0.1, num_of_exps, algo)

  print('H-simplex40')
  PolyList = GenSimplex(40, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-simplex40', (1/prod(1:40)), 0.1, num_of_exps, algo)

  print('H-simplex50')
  PolyList = GenSimplex(50, 'H')
  Hruntest(PolyList$A, PolyList$b, 'H-simplex50', (1/prod(1:50)), 0.1, num_of_exps, algo)

  if(algo=="SOB"){
    print('----------------------------------------')
    print('--------6th test [skinny_cubes]---------')
    print('----------------------------------------')
    cat('\n')
    print('H-skinny_cube10')
    PolyList = GenSkinnyCube(10)
    Hruntest(PolyList$A, PolyList$b, 'H-skinny_cube10', 102400, 0.1, num_of_exps, algo)
  
    print('H-skinny_cube20')
    PolyList = GenSkinnyCube(20)
    Hruntest(PolyList$A, PolyList$b, 'H-skinny_cube20', 104857600, 0.1, num_of_exps, algo)
  }
}
