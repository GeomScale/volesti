context("Volume Z test")

library(volesti)

Zruntest <- function(Mat, name_string, tol, num_of_exps, algo){
  
  exactvol = ExactZonoVol(Mat)
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "SOB") {
      vol = vol + volume(G=Mat, rounding=TRUE)
    } else {
      vol = vol + volume(G=Mat, CG=TRUE, error=0.2, rounding=TRUE)
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
    res = 1
  }
  test_that("Z volume test", {
    expect_equal(res, 1)
  })
}


num_of_exps = 10

for (i in 1:2) {
  if (i==1) {
    algo = 'CG'
    tol = 0.2
  } else {
    algo = 'SOB'
    tol = 0.1
  }

  ZonoMat = GenZonotope(2, 4)
  Zruntest(ZonoMat, 'Zonotope_2_4', tol, num_of_exps, algo)
  
  ZonoMat = GenZonotope(2, 8)
  Zruntest(ZonoMat, 'Zonotope_2_8', tol, num_of_exps, algo)
  

  
  ZonoMat = GenZonotope(4, 8)
  Zruntest(ZonoMat, 'Zonotope_4_8', tol, num_of_exps, algo)
  
  ZonoMat = GenZonotope(4, 10)
  Zruntest(ZonoMat, 'Zonotope_4_10', tol, num_of_exps, algo)
  
  
  
  ZonoMat = GenZonotope(5, 10)
  Zruntest(ZonoMat, 'Zonotope_5_10', tol, num_of_exps, algo)
  
}
