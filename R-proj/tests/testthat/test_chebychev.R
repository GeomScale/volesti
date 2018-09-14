context("Chebychev test")

library(volesti)


runCheTest <- function(A, b, name_string, radius, tol) {
  dimension = dim(A)[2]
  vec_ball = CheBall(A,b)
  rad = vec_ball[dimension + 1]
  
  error = abs(radius - rad) / radius
  #print(paste0('error = ',error))
  if (error >= tol){
    #print(paste0('TEST FAILED!! ', error, ' > ', tol))
    res = 0
  } else {
    #print(paste0('Test PASSED!! [', name_string, ']'))
    res = 1
  }
  print(res)
  #cat('\n')
  return(res)
}


path = system.file('extdata', package = 'volesti')
tol = 0.00001

PolyList = GenCube(10, 'H')
runCheTest(PolyList$A, PolyList$b, 'H-cube10', 1.0, tol)

PolyList = GenCube(20, 'H')
test_that("Chebychev test", {
  res = runCheTest(PolyList$A, PolyList$b, 'H-cube20', 1.0, tol)
  expect_equal(res, 1)
})


PolyList = GenCube(30, 'H')

test_that("Chebychev test", {
  res = runCheTest(PolyList$A, PolyList$b, 'H-cube30', 1.0, tol)
  expect_equal(res, 1)
})



PolyList = GenCross(10, 'H')
res = runCheTest(PolyList$A, PolyList$b, 'H-cross10', 0.316228, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})



ListPoly = fileToMatrix(paste0(path,'/birk3.ine'))
res = runCheTest(ListPoly$A, ListPoly$b, 'H-birk3', 0.207107, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

ListPoly = fileToMatrix(paste0(path,'/birk4.ine'))
res = runCheTest(ListPoly$A, ListPoly$b, 'H-birk4', 0.122008, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

ListPoly = fileToMatrix(paste0(path,'/birk5.ine'))
res = runCheTest(ListPoly$A, ListPoly$b, 'H-birk5', 0.0833333, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

ListPoly = fileToMatrix(paste0(path,'/birk6.ine'))
res = runCheTest(ListPoly$A, ListPoly$b, 'H-birk6', 0.0618034, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})


PolyList = GenProdSimplex(5)
res = runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_5_5', 0.138197, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenProdSimplex(10)
res = runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_10_10', 0.0759747, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenProdSimplex(15)
res = runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_15_15', 0.0529858, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenProdSimplex(20)
res = runCheTest(PolyList$A, PolyList$b, 'H-prod_simplex_20_20', 0.0408628, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})


PolyList = GenSimplex(10, 'H')
res = runCheTest(PolyList$A, PolyList$b, 'H-simplex10', 0.0759747, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenSimplex(20, 'H')
res = runCheTest(PolyList$A, PolyList$b, 'H-simplex20', 0.0408628, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenSimplex(30, 'H')
res = runCheTest(PolyList$A, PolyList$b, 'H-simplex30', 0.0281871, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenSimplex(40, 'H')
res = runCheTest(PolyList$A, PolyList$b, 'H-simplex40', 0.0215868, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenSimplex(50, 'H')
res = runCheTest(PolyList$A, PolyList$b, 'H-simplex50', 0.017522, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})


PolyList = GenSkinnyCube(10)
res = runCheTest(PolyList$A, PolyList$b, 'H-skinny_cube10', 1.0, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})

PolyList = GenSkinnyCube(20)
res = runCheTest(PolyList$A, PolyList$b, 'H-skinny_cube20', 1.0, tol)
test_that("Chebychev test", {
  expect_equal(res, 1)
})
