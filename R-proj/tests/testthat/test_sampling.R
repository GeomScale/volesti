context("Sampling test")

library(volesti)

runsample <- function(P, name_string, dist){
  if (dist == "uniform") {
    p = sample_points(P)
  } else {
    p = sample_points(P, distribution = "gaussian")
  }
  if (length(p[is.nan(p)])>0 | length(p[is.infinite(p)])>0) {
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
    distribution = 'gaussian'
  } else {
    distribution = 'uniform'
  }
  
  test_that("Sampling test", {
    P= GenCube(10, 'H')
    res = runsample(P, 'H-cube10', distribution)
    expect_equal(res, 1)
  })
  
  if (!cran_only) {
  test_that("Sampling test", {
    skip_on_cran()
    P = GenCube(20, 'H')
    res = runsample(P, 'H-cube20', distribution)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Sampling test", {
    skip_on_cran()
    P = GenCube(5, 'V')
    res = runsample(P, 'V-cube5', distribution)
    expect_equal(res, 1)
  })
  }

  test_that("Sampling test", {
    P = GenCross(10, 'H')
    res = runsample(P, 'H-cross10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = fileToMatrix(paste0(path,'/birk3.ine'))
    res = runsample(P, 'H-birk3', distribution)
    expect_equal(res, 1)
  })#
  
  if (!cran_only) {
  test_that("Sampling test", {
    skip_on_cran()
    P = fileToMatrix(paste0(path,'/birk4.ine'))
    res = runsample(P, 'H-birk4', distribution)
    expect_equal(res, 1)
  })
  }
  
  test_that("Sampling test", {
    P = GenProdSimplex(5)
    res = runsample(P, 'H-prod_simplex_5_5', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = GenProdSimplex(10)
    res = runsample(P, 'H-prod_simplex_10_10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = GenSimplex(10, 'H')
    res = runsample(P, 'H-prod_simplex10', distribution)
    expect_equal(res, 1)
  })
  
  if (!cran_only) {
  test_that("Sampling test", {
    skip_on_cran()
    P = GenSimplex(20, 'H')
    res = runsample(P, 'H-simplex20', distribution)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Sampling test", {
    skip_on_cran()
    P = GenSimplex(10, 'V')
    res = runsample(P, 'V-simplex10', distribution)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Sampling test", {
    skip_on_cran()
    P = GenSimplex(20, 'V')
    res = runsample(P, 'V-simplex20', distribution)
    expect_equal(res, 1)
  })
  }
  
  test_that("Sampling test", {
    P = GenSkinnyCube(10)
    res = runsample(P, 'H-skinny_cube10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    P = GenSkinnyCube(20)
    res = runsample(P, 'H-skinny_cube20', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    Z = GenZonotope(4, 8)
    res = runsample(Z, 'zonotope_4_8', distribution)
    expect_equal(res, 1)
  })
  
}
