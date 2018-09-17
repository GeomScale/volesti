context("Rounding test")

library(volesti)

testRound <- function(Mat, vector, exactvol, tol, name_string, num_of_exps, algo, rotation){
  
  if (rotation) {
    listHpoly = rand_rotate(A=Mat, b=vector)
    listHpoly = round_polytope(A=listHpoly$A, b=listHpoly$b)
  } else {
    listHpoly = round_polytope(A=Mat, b=vector)
  }
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + listHpoly$round_value * volume(A=listHpoly$A, b=listHpoly$b)
    } else {
      vol = vol + listHpoly$round_value * volume(A=listHpoly$A, b=listHpoly$b, CG=TRUE, error=0.1)
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

for (i in 1:2) {
  if (i==1) {
    algo = 'CG'
    num_of_exps = 10
  } else {
    algo = 'SOB'
    num_of_exps = 10
  }
  
  #test_that("Rounding test", {
    #PolyList = GenSkinnyCube(10)
    #res = testRound(PolyList$A, PolyList$b, 102400, 0.1, 'H-skinny_cube10', num_of_exps, algo, TRUE)
    #expect_equal(res, 1)
  #})
  
  #PolyList = GenSkinnyCube(20)
  #testRound(PolyList$A, PolyList$b, 104857600, 0.3, 'H-skinny_cube20', num_of_exps, algo, TRUE)
  
  test_that("Rounding H-skinny_cube10", {
    PolyList = GenSkinnyCube(10)
    res = testRound(PolyList$A, PolyList$b, 102400, 0.1, 'H-skinny_cube10', num_of_exps, algo, FALSE)
    expect_equal(res, 1)
  })
  
  if (!cran_only) {
    skip_on_cran()
    test_that("Rounding H-skinny_cube20", {
      PolyList = GenSkinnyCube(20)
      res = testRound(PolyList$A, PolyList$b, 104857600, 0.3, 'H-skinny_cube20', num_of_exps, algo, FALSE)
      expect_equal(res, 1)
    })
  }
  
}
