library(volesti)


r10=c()
r15=c()
for (i in c(5,10,15,50)) {
  
  Z1=GenZonotope(10,i*10)
  st=vol_hzono(Z1,pca_ratio = TRUE)
  
  Z2=GenZonotope(15,i*15)
  st2=vol_hzono(Z2,pca_ratio = TRUE)
  
  r10=c(r10,(st[5]/st[1])^(1/10))
  r15=c(r10,(st2[5]/st2[1])^(1/15))
  
  save(r10, file = "pca_10_10_15_50.RData")
  save(r15, file = "pca_15_10_15_50.RData")
  
}