
library(volesti)
library(TruncatedNormal)

dims = c(2, seq(from=5, to=50, by=5))
time_dims = c()
for (i in dims){
  
  Z=GenZonotope(i,2*i)
  #print(Z)
  sum_tim = 0
  for (j in 1:10){
    tim=proc.time()
    test_vol = vol_zono(P = Z, e = 0.1, mvrandn =  mvrandn, verbose = FALSE)
    tim=proc.time()-tim
    sum_tim = sum_tim + as.numeric(as.character(tim[3]))
  }
  sum_tim = sum_tim / 10
  
  time_dims = c(time_dims ,sum_tim)
  
}
  