library(volesti)
#library(ggplot2)
#library(tidyr)

bodies1=c()
#bodies2=c()
continue_ball = TRUE
time_limit = 1000
i=2
while(TRUE) {
  print(paste0('i = ',i))
  if (!continue_ball) {
    break
  }
  num_b=0
  #num_b2 = 0
  for (j in 1:1){
    Z = GenZonotope(i,floor(1.4*i))
    x=system.time({ b1 = vol_hzono(Z,steps_only = TRUE)})
    if (as.numeric(x[3])>time_limit) {
      continue_ball = TRUE
    }
    num_b = num_b + b1
    #if (i<21) {
    #x=system.time({ b1 = cg_volume(Z,steps_only = TRUE)})
    #num_b2 = num_b2 + b1
    #}
  }
  
  num_b = num_b / 1
  #num_b2 = num_b2 / 5
  #bodies2=c(bodies2, num_b2)
  bodies1=c(bodies1, num_b)
  save(bodies1, file = "Hbodies_d_14d.RData")
  #save(bodies2, file = "bodies_d_15d.RData")
  if(i==2){
    i=5
  } else {
    i = i + 5
  }
}

#bodies2 = c(bodies2, rep(bodies2[length(bodies2)],length(bodies1)-length(bodies2)))

#test2.data <- data.frame(
#  dimension= 2:(i-1),
#  BallAnnealing = bodies1, 
#  CoolingGaussian = bodies2
#)

#save(test2.data, file = "d_15d_data.RData")

#test2.data %>%
#  gather(methods,bodies,BallAnnealing, CoolingGaussian) %>%
#  ggplot(aes(x=dimension, y=bodies, colour=methods)) +
#  geom_line() + geom_point()

