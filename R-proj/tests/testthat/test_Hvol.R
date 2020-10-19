context("H-polytopes' volume test")

library(volesti)


Hruntest <- function(P, name_string, exactvol, tol, num_of_exps, alg, seed){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (alg == "CB") {
      vol = vol + volume(P, settings = list("algorithm" = "CB", "seed" = seed))
    } else if (alg == "SOB") {
      vol = vol + volume(P, settings = list("algorithm" = "SOB", "seed" = seed))
    } else {
      vol = vol + volume(P, settings = list("algorithm" = "CG", "seed" = seed))
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
path = system.file('extdata', package = 'volesti')

for (i in 1:2) {
  
  if (i==1) {
    algo = 'CG'
    num_of_exps = 10
  } else {
    algo = 'CB'
    num_of_exps = 10
  }

  
  test_that("Volume H-cube10", {
    P = gen_cube(10, 'H')
    res = Hruntest(P, 'H-cube10', 1024, 0.2, num_of_exps, algo, 5)
    expect_equal(res, 1)
  })
  
  test_that("Volume H-cross5", {
    P = gen_cross(5, 'H')
    res = Hruntest(P, 'H-cross10', 0.2666667, 0.2, num_of_exps, algo, 5)
    expect_equal(res, 1)
  })
  

  test_that("Volume H-prod_simplex_5_5", {
    P = gen_prod_simplex(5)
    res = Hruntest(P, 'H-prod_simplex_5_5', (1/prod(1:5))^2, 0.2, num_of_exps, algo, 5)
    expect_equal(res, 1)
  })
  
  test_that("Volume H-cube10", {
    P = gen_cube(10, 'H')
    res = Hruntest(P, 'H-cube10', 1024, 0.2, 5, "SOB", 5)
    expect_equal(res, 1)
  })
  
}
