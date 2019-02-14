context("Zonotopes' volume test")

library(volesti)

Zruntest <- function(P, name_string, tol, num_of_exps, algo){
  
  exactvol = exact_vol(P)
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + volume(P, rounding=TRUE)
    } else {
      vol = vol + volume(P, error=0.2, Algo = "CG", rounding=TRUE)
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
    tol = 0.4
  } else {
    algo = 'SOB'
    tol = 0.2
  }

  test_that("Volume Zonotope_2_4", {
    Z = GenZonotope(2, 4)
    res = Zruntest(Z, 'Zonotope_2_4', tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  
  if (!cran_only) {
  test_that("Volume Zonotope_2_8", {
    skip_on_cran()
    Z = GenZonotope(2, 8)
    res = Zruntest(Z, 'Zonotope_2_8', tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Volume Zonotope_4_8", {
    skip_on_cran()
    Z = GenZonotope(4, 8)
    res = Zruntest(Z, 'Zonotope_4_8', tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Volume Zonotope_4_10", {
    skip_on_cran()
    Z = GenZonotope(4, 10)
    res = Zruntest(Z, 'Zonotope_4_10', tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Volume Zonotope_5_10", {
    skip_on_cran()
    Z = GenZonotope(5, 10)
    res = Zruntest(Z, 'Zonotope_5_10', tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
}
