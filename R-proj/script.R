library(volesti)
library(TruncatedNormal)
library(tmvtnorm)

Z=GenZonotope(5,10)
test_vol = vol_zono(Z,0.1,SampleGibbs, mvrandn, mvNcdf, TRUE,FALSE)