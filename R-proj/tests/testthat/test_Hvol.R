context("H-Volume test")

library(volesti)



Hruntest <- function(Mat, vector, name_string, exactvol, tol, num_of_exps, alg){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (alg == "SOB") {
      vol = vol + volume(A=Mat, b=vector)
    } else {
      vol = vol + volume(A=Mat, b=vector, CG=TRUE, error=0.1)
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

  
  test_that("Volume test", {
    PolyList = GenCube(10, 'H')
    res = Hruntest(PolyList$A, PolyList$b, 'H-cube10', 1024, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenCube(20, 'H')
      res = Hruntest(PolyList$A, PolyList$b, 'H-cube20', 1048576, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenCube(30, 'H')
      res = Hruntest(PolyList$A, PolyList$b, 'H-cube30', 1073742000, 0.2, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }
  
  test_that("Volume test", {
    PolyList = GenCross(10, 'H')
    res = Hruntest(PolyList$A, PolyList$b, 'H-cross10', 0.0002821869, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })

  test_that("Volume test", {
    ListPoly = fileToMatrix(paste0(path,'/birk3.ine'))
    res = Hruntest(ListPoly$A, ListPoly$b, 'H-birk3', 0.125, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
  
  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      ListPoly = fileToMatrix(paste0(path,'/birk4.ine'))
      res = Hruntest(ListPoly$A, ListPoly$b, 'H-birk4', 0.000970018, 0.2, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      ListPoly = fileToMatrix(paste0(path,'/birk5.ine'))
      res = Hruntest(ListPoly$A, ListPoly$b, 'H-birk5', 0.000000225, 0.2, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      ListPoly = fileToMatrix(paste0(path,'/birk6.ine'))
      res = Hruntest(ListPoly$A, ListPoly$b, 'H-birk6', 0.0000000000009455459196, 0.5, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  test_that("Volume test", {
    PolyList = GenProdSimplex(5)
    res = Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_5_5', (1/prod(1:5))^2, 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenProdSimplex(10)
      res = Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_10_10', (1/prod(1:10))^2, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenProdSimplex(15)
      res = Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_15_15', (1/prod(1:15))^2, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenProdSimplex(20)
      res = Hruntest(PolyList$A, PolyList$b, 'H-prod_simplex_20_20', (1/prod(1:20))^2, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }
    
  test_that("Volume test", {
    PolyList = GenSimplex(10, 'H')
    res = Hruntest(PolyList$A, PolyList$b, 'H-simplex10', (1/prod(1:10)), 0.1, num_of_exps, algo)
    expect_equal(res, 1)
  })
    
  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenSimplex(20, 'H')
      res = Hruntest(PolyList$A, PolyList$b, 'H-simplex20', (1/prod(1:20)), 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenSimplex(30, 'H')
      res = Hruntest(PolyList$A, PolyList$b, 'H-simplex30', (1/prod(1:30)), 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenSimplex(40, 'H')
      res = Hruntest(PolyList$A, PolyList$b, 'H-simplex40', (1/prod(1:40)), 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume test", {
      PolyList = GenSimplex(50, 'H')
      res = Hruntest(PolyList$A, PolyList$b, 'H-simplex50', (1/prod(1:50)), 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if(algo=="SOB"){

    test_that("Volume test", {
      PolyList = GenSkinnyCube(10)
      res = Hruntest(PolyList$A, PolyList$b, 'H-skinny_cube10', 102400, 0.1, num_of_exps, algo)
      expect_equal(res, 1)
    })
  
    if (!cran_only) {
      skip_on_cran()
      test_that("Volume test", {
        PolyList = GenSkinnyCube(20)
        res = Hruntest(PolyList$A, PolyList$b, 'H-skinny_cube20', 104857600, 0.1, num_of_exps, algo)
        expect_equal(res, 1)
      })
    }
  }
}
