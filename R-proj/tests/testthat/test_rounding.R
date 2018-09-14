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
  #print(paste0('volume approximation of ', name_string,' = ',vol))
  #print(paste0('exact volume of ', name_string,' = ',exactvol))
  error = abs(vol - exactvol) / exactvol
  #print(paste0('error = ',error))
  if (error >= tol){
    res = 0
  } else {
    res = 1
  }
  
  test_that("Rounding test", {
    expect_equal(res, 1)
  })
}


for (i in 1:2) {
  if (i==1) {
    algo = 'CG'
    num_of_exps = 10
  } else {
    algo = 'SOB'
    num_of_exps = 20
  }
  
  PolyList = GenSkinnyCube(10)
  testRound(PolyList$A, PolyList$b, 102400, 0.1, 'H-skinny_cube10', num_of_exps, algo, TRUE)
  
  PolyList = GenSkinnyCube(20)
  testRound(PolyList$A, PolyList$b, 104857600, 0.3, 'H-skinny_cube20', num_of_exps, algo, TRUE)
  
  PolyList = GenSkinnyCube(10)
  testRound(PolyList$A, PolyList$b, 102400, 0.1, 'H-skinny_cube10', num_of_exps, algo, FALSE)
  
  PolyList = GenSkinnyCube(20)
  testRound(PolyList$A, PolyList$b, 104857600, 0.3, 'H-skinny_cube20', num_of_exps, algo, FALSE)
  
}
