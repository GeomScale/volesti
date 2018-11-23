library(volesti)

tim=c()
vols=c()

for (i in c(2,seq(from=5,to=50,by=5))){
  
  print(paste0("i = ",i))
  Z=GenZonotope(i,4*i)
  x=system.time({ b1 = ban_volume(Z)})
  
  tim = c(tim, as.numeric(x[3]))
  #print(tim)
  vols = c(vols, b1[1])
  
  save(tim, file = "times_2_50_d_4d.RData")
  save(vols, file = "vols_2_50_d_4d.RData")
  
}

