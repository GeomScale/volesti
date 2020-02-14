context("Zonotopes' volume test")

library(volesti)

Zruntest <- function(P, name_string, tol, num_of_exps, algo){
  
  exactvol = exact_vol(P)
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "CB") {
      vol = vol + volume(P, rounding=TRUE)
    } else {
      vol = vol + volume(P, error=0.1, algo = "CG", rounding=TRUE)
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
    algo = 'CB'
    tol = 0.3
  }

  test_that("Volume Zonotope_2_4", {
    Z = gen_rand_zonotope(2, 4)
    res = Zruntest(Z, 'Zonotope_2_4', tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  
}
