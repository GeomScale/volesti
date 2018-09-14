context("Volume V test")

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
  #print(paste0('volume approximation of ',name_string,' = ',vol))
  #print(paste0('exact volume of ',name_string,' = ',exactvol))
  error = abs(vol - exactvol) / exactvol
  #print(paste0('error = ',error))
  if (error >= tol){
    res = 0
  } else {
    res =1
  }
  test_that("V volume test", {
    expect_equal(res, 1)
  })
}


num_of_exps = 10

for (i in 1:2) {
  
  if (i==1) {
    algo = 'CG'
    tol = 0.2
    num_of_exps = 20
  } else {
    algo = 'SOB'
    tol = 0.1
  }
  

  PolyMat = GenCube(3, 'V')
  Vruntest(PolyMat, 'V-cube3', 8, tol, num_of_exps, algo)

  PolyMat = GenCube(4, 'V')
  Vruntest(PolyMat, 'V-cube4', 16, tol, num_of_exps, algo)
  

  PolyMat = GenCross(5, 'V')
  Vruntest(PolyMat, 'V-cross5', 0.2666667, tol, num_of_exps, algo)
  

  PolyMat = GenCross(7, 'V')
  Vruntest(PolyMat, 'V-cross7', 0.02539683, tol, num_of_exps, algo)
  

  PolyMat = GenSimplex(5, 'V')
  Vruntest(PolyMat, 'V-simplex5', 1/prod(1:5), tol, num_of_exps, algo)

  PolyMat = GenSimplex(7, 'V')
  Vruntest(PolyMat, 'V-simplex7', 1/prod(1:7), tol, num_of_exps, algo)

}
