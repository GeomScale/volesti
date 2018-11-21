library(volesti)
#library(ggplot2)
#library(tidyr)

Z=GenZonotope(100,200)
x=system.time({ b1 = ban_volume(Z,e=0.2)})
save(x, file = "tim_100_200.RData")
save(b1, file = "vec_vol_100_200.RData")
