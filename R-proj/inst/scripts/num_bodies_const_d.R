library(volesti)
library(ggplot2)
library(tidyr)

mat_bodies = matrix(0,nrow = 11, ncol = 25)
count_dim = 0
for (i in seq(from=10,to=60, by=5)) {
  num_b=0
  count_dim = count_dim + 1
  count_gen = 0
  for (j in seq(from=i+1,to=3*i+1,by=5)){ 
    Z=GenZonotope(i,j)
    b1 = ban_volume(Z,steps_only = TRUE)
    count_gen = count_gen + 1
    mat_bodies[count_dim, count_gen] = b1
    save(mat_bodies, file = "mat_bodies.RData")
  }
}