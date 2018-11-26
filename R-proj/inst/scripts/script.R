library(volesti)
library(ggplot2)
library(tidyr)

time1=c()
time2=c()
steps1=c()
steps2=c()
path = system.file('extdata', package = 'volesti')
for (i in 2:200) {
  print(i)
  #name_bir = paste0('/birk',i,'.ine')
  #HP = fileToMatrix(paste0(path,name_bir))
  HP = GenCube(i,'H')
  tim=proc.time()
  st1=ban_volume(HP)
  tim=proc.time()-tim
  tim = as.numeric(as.character(tim[3]))
  time1=c(time1, tim)
  steps1=c(steps1, st1[3])
  
  tim=proc.time()
  st2=cg_volume(HP)
  tim=proc.time()-tim
  tim = as.numeric(as.character(tim[3]))
  time2=c(time2, tim)
  steps2=c(steps2, st2[3])
}

test.data <- data.frame(
  dimension= 2:200,
  CoolingBall = time1, 
  CoolingGaussian = time2
)

test2.data <- data.frame(
  dimension= 2:200,
  CoolingBall = steps1, 
  CoolingGaussian = steps2
)
save(test.data, file = "times_cubes_2_200.RData")
save(test2.data, file = "steps_cubes_2_200.RData")

test.data %>%
  gather(methods,steps, BallAnnealing, CoolingGaussian) %>%
  ggplot(aes(x=dimension, y=steps, colour=methods)) +
  geom_line() + geom_point()