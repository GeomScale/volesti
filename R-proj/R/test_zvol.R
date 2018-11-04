test_zvol <- function(d,m,verbose=TRUE){
  
  Z=GenZonotope(d,m)
  #print(Z)
  
  tim=proc.time()
  test_vol = vol_zono(Z,0.05,SampleGibbs, mvrandn, mvNcdf, verbose,FALSE)
  tim=proc.time()-tim
  print(paste0('test method volume = ',test_vol))
  print(paste0('test method time: ',as.numeric(as.character(tim[3]))))
  print(' ')
  
  #tim=proc.time()
  #ev = exact_vol(Z)
  #tim=proc.time()-tim
  #print(paste0('exact volume = ',ev))
  #print(paste0('exact volume time: ',as.numeric(as.character(tim[3]))))
  #print(' ')
  #ev=1
  
  tim=proc.time()
  est_vol = volume(P=Z, error = 0.1,Algo = list("CG"=TRUE), rounding = TRUE)
  tim=proc.time()-tim
  print(paste0('cg-HnR volume = ',est_vol))
  print(paste0('cg-HnR volume time: ',as.numeric(as.character(tim[3]))))
  print(' ')
  
  #print(paste0('test method error = ',abs(test_vol-ev)/ev))
  #print(paste0('cg-HnR error = ',abs(est_vol-ev)/ev))
  
  return(Z)
  
}