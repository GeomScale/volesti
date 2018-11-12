library(volesti)
library(ggplot2)
library(tidyr)

bodies1=c()
bodies2=c()
continue_ball = TRUE
continue_cg = TRUE
end_dim_cg = TRUE
#end_dim_ban = TRUE
time_limit = 1800
i=2
while(TRUE) {
  if (!continue_cg && !continue_ball) {
    break
  }
  print(paste0(i," ",continue_ball," ",continue_cg, " ",end_dim_cg))
  num_b=0
  num_b2 = 0
  #Z = GenZonotope(i,2*i)
  for (j in 1:5){
    Z = GenZonotope(i,2*i)
    x=system.time({ b1 = ban_volume(Z,steps_only = TRUE)})
    if (as.numeric(x[3])>time_limit) {
      continue_ball = FALSE
    }
    num_b = num_b + b1
    if (continue_cg || end_dim_cg) {
      x=system.time({ b2 = cg_volume(Z,steps_only = TRUE)})
      if (as.numeric(x[3])>time_limit) {
        continue_cg = FALSE
      }
      num_b2 = num_b2 + b2
    } else {
      num_b2 = 0
    }
  }
  num_b = num_b / 5
  num_b2 = num_b2 / 5
  bodies2=c(bodies2, num_b2)
  bodies1=c(bodies1, num_b)
  i = i + 1
  if (!continue_cg) {
    end_dim_cg = FALSE
  }
}
#bodies2 = c(bodies2, rep(bodies2[length(bodies2)],length(bodies1)-length(bodies2)))

test2.data <- data.frame(
  dimension= 2:(i-1),
  BallAnnealing = bodies1, 
  CoolingGaussian = bodies2
)

save(test2.data, file = "script2_data.RData")


dim = 51
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



#test2.data %>%
#  gather(methods,bodies,BallAnnealing, CoolingGaussian) %>%
#  ggplot(aes(x=dimension, y=bodies, colour=methods)) +
#  geom_line() + geom_point()

