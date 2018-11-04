SampleGibbs <- function(N, mu, sig, up, lb, burns, W){
  
  X2 <- t(rtmvnorm(n=N, mean = mu, sigma=sig,lower=lb, upper=up, algorithm="gibbs", burn.in.samples=burns, start.value=mu, thinning=W))
  print(X2[,1])
  
  return(X2)
  
}