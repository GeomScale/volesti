
library(volesti)

d=100
m=100
P=GenCube(3,'H')

tim = system.time({ vol = volume(P, nn=d, mm=m, rounding = FALSE) })

print(vol)
print(tim)