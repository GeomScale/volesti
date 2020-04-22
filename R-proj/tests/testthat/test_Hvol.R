context("H-polytopes' volume test")

library(volesti)


Hruntest <- function(P, name_string, exactvol, tol, num_of_exps, alg){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (alg == "SOB") {
      vol = vol + volume(P)
    } else {
      vol = vol + volume(P, Algo = "CG")
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
    num_of_exps = 20
  } else {
    algo = 'SOB'
    num_of_exps = 10
  }

  
  test_that("Volume H-cube10", {
    P = GenCube(10, 'H')
    res = Hruntest(P, 'H-cube10', 1024, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
    test_that("Volume H-cube20", {
      skip_on_cran()
      P = GenCube(20, 'H')
      res = Hruntest(P, 'H-cube20', 1048576, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

      
  if (!cran_only) {
  test_that("Volume H-cube30", {
    skip_on_cran()
    P = GenCube(30, 'H')
    res = Hruntest(P, 'H-cube30', 1073742000, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  test_that("Volume H-cross5", {
    P = GenCross(5, 'H')
    res = Hruntest(P, 'H-cross10', 0.2666667, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
  test_that("Volume H-birk3", {
    skip_on_cran()
    P = fileToMatrix(paste0(path,'/birk3.ine'))
    res = Hruntest(P, 'H-birk3', 0.125, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Volume H-birk4", {
    skip_on_cran()
    P = fileToMatrix(paste0(path,'/birk4.ine'))
    res = Hruntest(P, 'H-birk4', 0.000970018, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-birk5", {
    skip_on_cran()
    P = fileToMatrix(paste0(path,'/birk5.ine'))
    res = Hruntest(P, 'H-birk5', 0.000000225, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-birk6", {
    skip_on_cran()
    P = fileToMatrix(paste0(path,'/birk6.ine'))
    res = Hruntest(P, 'H-birk6', 0.0000000000009455459196, 0.5, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  test_that("Volume H-prod_simplex_5_5", {
    P = GenProdSimplex(5)
    res = Hruntest(P, 'H-prod_simplex_5_5', (1/prod(1:5))^2, 0.2, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
  test_that("Volume H-prod_simplex_10_10", {
    skip_on_cran()
    P = GenProdSimplex(10)
    res = Hruntest(P, 'H-prod_simplex_10_10', (1/prod(1:10))^2, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-prod_simplex_15_15", {
    skip_on_cran()
    P = GenProdSimplex(15)
    res = Hruntest(P, 'H-prod_simplex_15_15', (1/prod(1:15))^2, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-prod_simplex_20_20", {
    skip_on_cran()
    P = GenProdSimplex(20)
    res = Hruntest(P, 'H-prod_simplex_20_20', (1/prod(1:20))^2, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
    
  if (!cran_only) {
  test_that("Volume H-simplex10", {
    skip_on_cran()
    P = GenSimplex(10, 'H')
    res = Hruntest(P, 'H-simplex10', (1/prod(1:10)), 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }
  
  if (!cran_only) {
  test_that("Volume H-simplex20", {
    skip_on_cran()
    P = GenSimplex(20, 'H')
    res = Hruntest(P, 'H-simplex20', (1/prod(1:20)), 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-simplex30", {
    skip_on_cran()
    P = GenSimplex(30, 'H')
    res = Hruntest(P, 'H-simplex30', (1/prod(1:30)), 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-simplex40", {
    skip_on_cran()
    P = GenSimplex(40, 'H')
    res = Hruntest(P, 'H-simplex40', (1/prod(1:40)), 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if (!cran_only) {
  test_that("Volume H-simplex50", {
    skip_on_cran()
    P = GenSimplex(50, 'H')
    res = Hruntest(P, 'H-simplex50', (1/prod(1:50)), 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  }

  if(algo=="SOB"){

    if (!cran_only) {
    test_that("Volume H-skinny_cube10", {
      skip_on_cran()
      P = GenSkinnyCube(10)
      res = Hruntest(P, 'H-skinny_cube10', 102400, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
    }
  
    if (!cran_only) {
    test_that("Volume H-skinny_cube20", {
      skip_on_cran()
      P = GenSkinnyCube(20)
      res = Hruntest(P, 'H-skinny_cube20', 104857600, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
    }
    
  }
}
