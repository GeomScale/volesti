library(volesti)
library(TruncatedNormal)

gens = seq(from=15, to=90, by=5)

errors=t(matrix(0,6))
r = c(0,0,0,0,0,0)
for (i in gens){
  
  Z=GenZonotope(10,i)
  #print(Z)
  cg_vol = 0
  nm_vol = 0
  for (j in 1:5){
    test_vol = vol_zono(P = Z, e = 0.05, mvrandn =  mvrandn, verbose = FALSE)
    est_vol = volume(P=Z, error = 0.1,Algo = list("CG"=TRUE))
    cg_vol = cg_vol + est_vol
    nm_vol = nm_vol + test_vol
  }
  cg_vol = cg_vol / 5
  nm_vol = nm_vol / 5
  r[1] = cg_vol
  r[2] = cg_vol / 1.1
  r[3] = cg_vol / 0.9
  r[4] = nm_vol
  
  if (r[4]>r[2] && r[4]<r[3]) {
    r[5]=1
  } else {
    r[5] = 0
  }
  
  if (r[4]>(cg_vol / 1.2) && r[4]<(cg_vol / 0.8)) {
    r[6]=1
  } else {
    r[6] = 0
  }
  
  errors=rbind(errors, r)
  
}