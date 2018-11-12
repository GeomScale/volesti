library(volesti)
library(ggplot2)
library(tidyr)

dim =51
mat_data=matrix(0, nrow=dim-1, ncol=100)

for (i in 2:dim) {
  k=i+1
  m = min(100,4*k)
  for (j in k:m){
    num_b = 0
    for (l in 1:5) {
      Z=GenZonotope(i,j)
      num_b = num_b + ban_volume(Z,steps_only = TRUE)
    }
    num_b = num_b / 5
    mat_data[i-1,j] = num_b
  }
}
save(mat_data, file = "script3_data.RData")
