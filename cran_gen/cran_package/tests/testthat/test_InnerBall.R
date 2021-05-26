context("Chebychev ball test")

library(volesti)

runCheTest <- function(P, name_string, radius, tol) {
  
  vec_ball = inner_ball(P, lpsolve = TRUE)
  rad = vec_ball[length(vec_ball)]
  
  error = abs(radius - rad) / radius
  if (error >= tol){
    res = 0
  } else {
    res = 1
  }
  return(res)
}

tol = 0.00001

test_that("Chebychev test", {
  P = gen_cube(10, 'H')
  res = runCheTest(P, 'H-cube10', 1.0, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_cube(20, 'H')
  res = runCheTest(P, 'H-cube20', 1.0, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_cube(30, 'H')
  res = runCheTest(P, 'H-cube30', 1.0, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_cross(10, 'H')
  res = runCheTest(P, 'H-cross10', 0.316228, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_birkhoff(3)
  res = runCheTest(P, 'H-birk3', 0.207107, tol)
  expect_equal(res, 1)
})


test_that("Chebychev test", {
  P = gen_birkhoff(4)
  res = runCheTest(P, 'H-birk4', 0.122008, tol)
  expect_equal(res, 1)
})


test_that("Chebychev test", {
  P = gen_birkhoff(5)
  res = runCheTest(P, 'H-birk5', 0.0833333, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_birkhoff(6)
  res = runCheTest(P, 'H-birk6', 0.0618034, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_prod_simplex(5)
  res = runCheTest(P, 'H-prod_simplex_5_5', 0.138197, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_prod_simplex(10)
  res = runCheTest(P, 'H-prod_simplex_10_10', 0.0759747, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_prod_simplex(15)
  res = runCheTest(P, 'H-prod_simplex_15_15', 0.0529858, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_prod_simplex(20)
  res = runCheTest(P, 'H-prod_simplex_20_20', 0.0408628, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_simplex(10, 'H')
  res = runCheTest(P, 'H-simplex10', 0.0759747, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_simplex(20, 'H')
  res = runCheTest(P, 'H-simplex20', 0.0408628, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_simplex(30, 'H')
  res = runCheTest(P, 'H-simplex30', 0.0281871, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_simplex(40, 'H')
  res = runCheTest(P, 'H-simplex40', 0.0215868, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_simplex(50, 'H')
  res = runCheTest(P, 'H-simplex50', 0.017522, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_skinny_cube(10)
  res = runCheTest(P, 'H-skinny_cube10', 1.0, tol)
  expect_equal(res, 1)
})

test_that("Chebychev test", {
  P = gen_skinny_cube(20)
  res = runCheTest(P, 'H-skinny_cube20', 1.0, tol)
  expect_equal(res, 1)
})
