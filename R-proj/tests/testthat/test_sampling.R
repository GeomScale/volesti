context("Sampling test")

library(volesti)

runsample <- function(Mat, vector, Vpoly, Zono, name_string, dist){
  if (dist == "uniform") {
    if (Zono) {
      p = sample_points(G = Mat)
    } else if (Vpoly) {
      p = sample_points(V = Mat)
    } else {
      p = sample_points(A=Mat, b=vector)
    }
  } else {
    if (Zono) {
      p = sample_points(G=Mat, gaussian=TRUE)
    } else if (Vpoly) {
      p = sample_points(V=Mat, gaussian=TRUE)
    } else {
      p = sample_points(A=Mat, b=vector, gaussian=TRUE)
    }
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
    PolyList = GenCube(10, 'H')
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-cube10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenCube(20, 'H')
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-cube20', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyMat = GenCube(5, 'V')
    res = runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cube5', distribution)
    expect_equal(res, 1)
  })

  test_that("Sampling test", {
    PolyList = GenCross(10, 'H')
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-cross10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    ListPoly = fileToMatrix(paste0(path,'/birk3.ine'))
    res = runsample(ListPoly$A, ListPoly$b, FALSE, FALSE, 'H-birk3', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    ListPoly = fileToMatrix(paste0(path,'/birk4.ine'))
    res = runsample(ListPoly$A, ListPoly$b, FALSE, FALSE, 'H-birk4', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenProdSimplex(5)
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex_5_5', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenProdSimplex(10)
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex_10_10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenSimplex(10, 'H')
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenSimplex(20, 'H')
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-simplex20', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyMat = GenSimplex(10, 'V')
    res = runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-simplex10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyMat = GenSimplex(20, 'V')
    res = runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-simplex20', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenSkinnyCube(10)
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-skinny_cube10', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    PolyList = GenSkinnyCube(20)
    res = runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-skinny_cube20', distribution)
    expect_equal(res, 1)
  })
  
  test_that("Sampling test", {
    zonotope = GenZonotope(4, 8)
    res = runsample(zonotope, c(0,0), FALSE, TRUE, 'zonotope_4_8', distribution)
    expect_equal(res, 1)
  })
  
}
