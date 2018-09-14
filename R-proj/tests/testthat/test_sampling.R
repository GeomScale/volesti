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
  
  test_that("Sampling test", {
    expect_equal(res, 1)
  })
  
}


path = system.file('extdata', package = 'volesti')

for (i in 1:2) {
  
  if (i==1) {
    distribution = 'gaussian'
  } else {
    distribution = 'uniform'
  }
  

  PolyList = GenCube(10, 'H')
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-cube10', distribution)
  
  PolyList = GenCube(20, 'H')
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-cube20', distribution)
  
  PolyMat = GenCube(5, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cube5', distribution)
  
  PolyMat = GenCube(10, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cube10', distribution)

  PolyList = GenCross(10, 'H')
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-cross10', distribution)
  
  PolyMat = GenCross(20, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-cross20', distribution)
  

  ListPoly = fileToMatrix(paste0(path,'/birk3.ine'))
  runsample(ListPoly$A, ListPoly$b, FALSE, FALSE, 'H-birk3', distribution)
  
  ListPoly = fileToMatrix(paste0(path,'/birk4.ine'))
  runsample(ListPoly$A, ListPoly$b, FALSE, FALSE, 'H-birk4', distribution)
  
  ListPoly = fileToMatrix(paste0(path,'/birk5.ine'))
  runsample(ListPoly$A, ListPoly$b, FALSE, FALSE, 'H-birk5', distribution)
  
  ListPoly = fileToMatrix(paste0(path,'/birk6.ine'))
  runsample(ListPoly$A, ListPoly$b, FALSE, FALSE, 'H-birk6', distribution)
  

  PolyList = GenProdSimplex(5)
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex_5_5', distribution)
  
  PolyList = GenProdSimplex(10)
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex_10_10', distribution)
  
  PolyList = GenProdSimplex(15)
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex_15_15', distribution)
  
  PolyList = GenProdSimplex(20)
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex_20_20', distribution)
  

  PolyList = GenSimplex(10, 'H')
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-prod_simplex10', distribution)
  
  PolyList = GenSimplex(20, 'H')
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-simplex20', distribution)
  
  PolyMat = GenSimplex(10, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-simplex10', distribution)
  
  PolyMat = GenSimplex(20, 'V')
  runsample(PolyMat, c(0,0), TRUE, FALSE, 'V-simplex20', distribution)
  

  PolyList = GenSkinnyCube(10)
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-skinny_cube10', distribution)
  
  PolyList = GenSkinnyCube(20)
  runsample(PolyList$A, PolyList$b, FALSE, FALSE, 'H-skinny_cube20', distribution)
  

  zonotope = GenZonotope(4, 8)
  runsample(zonotope, c(0,0), FALSE, TRUE, 'zonotope_4_8', distribution)
  
  zonotope = GenZonotope(5, 10)
  runsample(zonotope, c(0,0), FALSE, TRUE, 'zonotope_5_10', distribution)
  
}
