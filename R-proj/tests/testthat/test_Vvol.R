context("V-polytopes' volume test")

library(volesti)

Vruntest <- function(P, name_string, exactvol, tol, num_of_exps, algorithm,seed){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (algorithm == "CB") {
      vol = vol + volume(P, rounding = "none", seed = seed)$volume
    } else {
      vol = vol + volume(P, settings = list("algorithm" = "CG", "error" = 0.1), rounding = "none", seed = seed)$volume
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
num_of_exps = 2

for (i in 1:2) {
  seed = 5
  if (i==1) {
    algo = 'CG'
    tol = 0.2
  } else {
    algo = 'CB'
    tol = 0.2
  }
  
  test_that("Volume V-simplex3", {
    P = gen_simplex(3, 'V')
    res = Vruntest(P, 'V-simplex3', 1/prod(1:3), tol, num_of_exps, algo, seed)
    expect_equal(res, 1)
  })

}
