context("Rounding test")

library(volesti)

testRound <- function(P, exactvol, tol, name_string, num_of_exps, algo, rotation){
  
  if (rotation) {
    P = rand_rotate(P)
    listHpoly = round_polytope(P)
  } else {
    listHpoly = round_polytope(P)
  }
  vol = 0
  for (j in 1:num_of_exps) {
    if (algo == "CB") {
      vol = vol + listHpoly$round_value * volume(listHpoly$P)
    } else {
      vol = vol + listHpoly$round_value * volume(listHpoly$P, error=0.1, algo=list("CG"=TRUE))
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
  
  num_of_exps = 10
  
  
  
  test_that("Rounding H-skinny_cube10", {
    P = gen_skinny_cube(10)
    res = testRound(P, 102400, 0.3, 'H-skinny_cube10', num_of_exps, 'CB', FALSE)
    expect_equal(res, 1)
  })
  

}
