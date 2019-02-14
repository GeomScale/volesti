context("V-polytopes' volume test")

library(volesti)

Vruntest <- function(P, name_string, exactvol, tol, num_of_exps, algo){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + volume(P)
    } else {
      vol = vol + volume(P, error=0.2, Algo = "CG")
    }
  }
  vol = vol / num_of_exps
  error = abs(vol - exactvol) / exactvol
  if (error >= tol){
    res = 0
  } else {
    res = 1
  }
  return(res)
}

cran_only = TRUE
num_of_exps = 5

for (i in 1:2) {
  
  if (i==1) {
    algo = 'CG'
    tol = 0.3
  } else {
    algo = 'SOB'
    tol = 0.2
  }
  
  if (!cran_only) {
  test_that("Volume V-cube3", {
    skip_on_cran()
    P = GenCube(3, 'V')
    res = Vruntest(P, 'V-cube3', 8, tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume V-cube4", {
    skip_on_cran()
    P = GenCube(4, 'V')
    res = Vruntest(P, 'V-cube4', 16, tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Volume V-cross3", {
    skip_on_cran()
    P = GenCross(3, 'V')
    res = Vruntest(P, 'V-cross3', 1.333333, tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("V volume test", {
    skip_on_cran()
    P = GenCross(7, 'V')
    Vruntest(P, 'V-cross7', 0.02539683, tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  test_that("Volume V-simplex3", {
    P = GenSimplex(3, 'V')
    res = Vruntest(P, 'V-simplex3', 1/prod(1:3), tol, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
  test_that("Volume V-simplex7", {
    skip_on_cran()
    P = GenSimplex(7, 'V')
    res = Vruntest(P, 'V-simplex7', 1/prod(1:7), tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

}
