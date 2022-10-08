context("Sampling test")

library(volesti)

runsample <- function(P, name_string, dist){
  if (dist == "uniform") {
    p = sample_points(P, n = 100)
  } else {
    p = sample_points(P, n = 100, distribution = list("density" = "gaussian"))
  }
  if (length(p[is.nan(p)])>0 | length(p[is.infinite(p)])>0) {
    res = 0
  } else {
    res = 1
  }
  return(res)
  
}

logconcave_sample<- function(P,distribution, n_samples ,n_burns){
if(distribution=="uniform"){
f <- function(x) (0)
grad_f <- function(x) (0)
L=1
m=1
pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f, "L_" = L, "m" = m))
return(psrf_univariate(pts))
}
else if(distribution== "gaussian"){
pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "variance"=8))
return(psrf_univariate(pts))
}
}

for (i in 1:2) {
  
  if (i==1) {
    distribution = 'gaussian'
  } else {
    distribution = 'uniform'
  }
  
  test_that("Sampling test", {
    P= gen_cube(10, 'H')
    res = runsample(P, 'H-cube10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = gen_cross(10, 'H')
    res = runsample(P, 'H-cross10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = gen_prod_simplex(5)
    res = runsample(P, 'H-prod_simplex_5_5', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = gen_prod_simplex(10)
    res = runsample(P, 'H-prod_simplex_10_10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = gen_simplex(10, 'H')
    res = runsample(P, 'H-prod_simplex10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = gen_skinny_cube(10)
    res = runsample(P, 'H-skinny_cube10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = gen_skinny_cube(20)
    res = runsample(P, 'H-skinny_cube20', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    Z = gen_rand_zonotope(4, 8)
    res = runsample(Z, 'zonotope_4_8', distribution)
    expect_equal(res, 1)
  })
  

  test_that("Sampling test", {
    P = gen_simplex(10, 'H')
    psrf = logconcave_sample(P,distribution,5000,2000)
    expect_lte(max(psrf), 1.1)
  })

}
