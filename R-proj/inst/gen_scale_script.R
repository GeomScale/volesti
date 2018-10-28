library(volesti)
library(TruncatedNormal)

gens = seq(from=15, to=90, by=5)

time_gens = c()
for (i in gens){
  
  Z=GenZonotope(10,i)
  #print(Z)
  sum_tim = 0
  for (j in 1:10){
    tim=proc.time()
    test_vol = vol_zono(P = Z, e = 0.05, mvrandn =  mvrandn, verbose = FALSE)
    tim=proc.time()-tim
    sum_tim = sum_tim + as.numeric(as.character(tim[3]))
  }
  sum_tim = sum_tim / 10
  
  time_gens = c(time_dims ,sum_tim)
  
}