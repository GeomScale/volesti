compute_indicators <- function(MatReturns, win_len = NULL, numSlices = NULL, N = NULL) {
  
  if (is.null(win_len)) win_len = 60
  if (is.null(numSlices)) numSlices = 100
  if (is.null(N)) N = 500000
  
  nrows = dim(MatReturns)[1]
  nassets = dim(MatReturns)[2]
  wl = win_len-1
  
  indicators = c()
  print(nrows-wl)
  for (i in 1:(nrows-wl)) {
    
    W=i:(i+wl)
    #nR = MatReturns[W,]
    E = cov(MatReturns[W,])
    
    compRet = rep(1,nassets)
    for (j in 1:nassets) {
      for (k in W) {
        compRet[j] = compRet[j] * (1 + MatReturns[k, j])
      }
      compRet[j] = compRet[j] - 1
    }
    
    cop = copula2(h = compRet, E = E, numSlices = numSlices, N = N)
    blue_mass = 0
    red_mass = 0
    
    for (row in 1:100) {
      for (col in 1:100) {
        if (row-col<=0.2*numSlices && row-col>=-0.2*numSlices) {
          if (row+col<0.8*numSlices || row+col>1.2*numSlices) {
            red_mass = red_mass + cop[row,col]
          }
        } else {
          if (row+col>=0.8*numSlices && row+col<=1.2*numSlices) {
            blue_mass = blue_mass + cop[row,col]
          }
        }
      }
    }
    indicators = c(indicators, blue_mass / red_mass)
    print(length(indicators))
  }
  
  return(indicators)
  
}