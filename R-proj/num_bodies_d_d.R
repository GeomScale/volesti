library(volesti)
library(ggplot2)
library(tidyr)

bodies1=c()
bodies2=c()
continue_ball = TRUE
time_limit = 1800
i=2
while(TRUE) {
  if (!continue_ball) {
    break
  }
  num_b=0
  num_b2 = 0
  for (j in 1:1){
    Z = GenZonotope(i,i+1)
    x=system.time({ b1 = ban_volume(Z,steps_only = TRUE)})
    if (as.numeric(x[3])>time_limit) {
      continue_ball = FALSE
    }
    num_b = num_b + b1
    if (i<41) {
      x=system.time({ b1 = cg_volume(Z,steps_only = TRUE)})
      num_b2 = num_b2 + b1
    }
  }
  
  num_b = num_b / 1
  num_b2 = num_b2 / 1
  bodies2=c(bodies2, num_b2)
  bodies1=c(bodies1, num_b)
  save(bodies1, file = "bodies1_d_d.RData")
  save(bodies2, file = "bodies2_d_d.RData")
  if(i==2){
    i=5
  } else {
    i = i + 5
  }
}

#bodies2 = c(bodies2, rep(bodies2[length(bodies2)],length(bodies1)-length(bodies2)))

#test2.data <- data.frame(
  #dimension= 2:(i-1),
  #BallAnnealing = bodies1, 
  #CoolingGaussian = bodies2
#)

#save(test2.data, file = "d_d_data.RData")

#test2.data %>%
#  gather(methods,bodies,BallAnnealing, CoolingGaussian) %>%
#  ggplot(aes(x=dimension, y=bodies, colour=methods)) +
#  geom_line() + geom_point()

