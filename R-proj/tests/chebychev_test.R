runCheTest <- function(A, b, name_string, radius, tol) {
  dimension = dim(A)[2]
  vec_ball = CheBall(A,b)
  rad = vec_ball[dimension + 1]
  
  error = abs(radius - rad) / radius
  print(paste0('error = ',error))
  if (error >= tol){
    print(paste0('TEST FAILED!! ', error, ' > ', tol))
  } else {
    print(paste0('Test PASSED!! [', name_string, ']'))
  }
  cat('\n')
}

library(volesti)

path = system.file('extdata', package = 'volesti')
tol = 0.00001

print('----------------------------------------')
print('------------1st test [cubes]------------')
print('----------------------------------------')
cat('\n')
print('H-cube10')
PolyList = GenCube(10, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-cube10', 1.0, tol)

print('H-cube20')
PolyList = GenCube(20, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-cube20', 1.0, tol)

print('H-cube30')
PolyList = GenCube(30, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-cube30', 1.0, tol)

print('----------------------------------------')
print('------2nd test [cross_polytopes]--------')
print('----------------------------------------')
cat('\n')

print('H-cross10')
PolyList = GenCross(10, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-cross10', 0.316228, tol)

print('----------------------------------------')
print('----------3rd test [birkhoff]-----------')
print('----------------------------------------')
cat('\n')
print('H-birk3')
ListPoly = fileToMatrix(paste0(path,'/birk3.ine'))
runCheTest(ListPoly$A, ListPoly$b, 'H-birk3', 0.207107, tol)

print('H-birk4')
ListPoly = fileToMatrix(paste0(path,'/birk4.ine'))
runCheTest(ListPoly$A, ListPoly$b, 'H-birk4', 0.122008, tol)

print('H-birk5')
ListPoly = fileToMatrix(paste0(path,'/birk5.ine'))
runCheTest(ListPoly$A, ListPoly$b, 'H-birk5', 0.0833333, tol)

print('H-birk6')
ListPoly = fileToMatrix(paste0(path,'/birk6.ine'))
runCheTest(ListPoly$A, ListPoly$b, 'H-birk6', 0.0618034, tol)

print('----------------------------------------')
print('--------4th test [prod_simplex]---------')
print('----------------------------------------')
cat('\n')
print('H-prod_simplex_5_5')
PolyList = GenProdSimplex(5)
runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_5_5', 0.138197, tol)

print('H-prod_simplex_10_10')
PolyList = GenProdSimplex(10)
runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_10_10', 0.0759747, tol)

print('H-prod_simplex_15_15')
PolyList = GenProdSimplex(15)
runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_15_15', 0.0529858, tol)

print('H-prod_simplex_20_20')
PolyList = GenProdSimplex(20)
runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_20_20', 0.0408628, tol)

print('----------------------------------------')
print('--------5th test [simplex]---------')
print('----------------------------------------')
cat('\n')
print('H-simplex10')
PolyList = GenSimplex(10, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-simplex10', 0.0759747, tol)

print('H-simplex20')
PolyList = GenSimplex(20, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-simplex20', 0.0408628, tol)

print('H-simplex30')
PolyList = GenSimplex(30, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-simplex30', 0.0281871, tol)

print('H-simplex40')
PolyList = GenSimplex(40, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-simplex40', 0.0215868, tol)

print('H-simplex50')
PolyList = GenSimplex(50, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-simplex50', 0.017522, tol)

print('----------------------------------------')
print('--------6th test [skinny_cubes]---------')
print('----------------------------------------')
cat('\n')
print('H-skinny_cube10')
PolyList = GenSkinnyCube(10)
runCheTest(PolyList$A, PolyList$b, 'H-skinny_cube10', 1.0, tol)

print('H-skinny_cube20')
PolyList = GenSkinnyCube(20)
runCheTest(PolyList$A, PolyList$b, 'H-skinny_cube20', 1.0, tol)
