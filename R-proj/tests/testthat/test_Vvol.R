context("V-polytopes' volume test")

library(volesti)

Vruntest <- function(Mat, name_string, exactvol, tol, num_of_exps, algo){
  
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + volume(V=Mat)
    } else {
      vol = vol + volume(V=Mat, CG=TRUE, error=0.2)
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
num_of_exps = 10

for (i in 1:2) {
  
  if (i==1) {
    algo = 'CG'
    tol = 0.2
  } else {
    algo = 'SOB'
    tol = 0.1
  }
  
  if (!cran_only) {
    skip_on_cran()
    test_that("Volume V-cube3", {
      PolyMat = GenCube(3, 'V')
      res = Vruntest(PolyMat, 'V-cube3', 8, tol, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume V-cube4", {
      PolyMat = GenCube(4, 'V')
      res = Vruntest(PolyMat, 'V-cube4', 16, tol, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }
  
  test_that("Volume V-cross3", {
    PolyMat = GenCross(3, 'V')
    res = Vruntest(PolyMat, 'V-cross3', 1.333333, tol, num_of_exps, algo)
    expect_equal(res, 1)
  })
  
  if (!cran_only) {
    skip_on_cran()
    test_that("V volume test", {
      PolyMat = GenCross(7, 'V')
      Vruntest(PolyMat, 'V-cross7', 0.02539683, tol, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }
  
  test_that("Volume V-simplex3", {
    PolyMat = GenSimplex(3, 'V')
    res = Vruntest(PolyMat, 'V-simplex3', 1/prod(1:3), tol, num_of_exps, algo)
    expect_equal(res, 1)
  })

  if (!cran_only) {
    skip_on_cran()
    test_that("Volume V-simplex7", {
      PolyMat = GenSimplex(7, 'V')
      res = Vruntest(PolyMat, 'V-simplex7', 1/prod(1:7), tol, num_of_exps, algo)
      expect_equal(res, 1)
    })
  }

}
