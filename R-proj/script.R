library(volesti)
library(ggplot2)
library(tidyr)

time1=c()
time2=c()
steps1=c()
steps2=c()
path = system.file('extdata', package = 'volesti')
for (i in 2:100) {
  print(i)
  #name_bir = paste0('/birk',i,'.ine')
  #HP = fileToMatrix(paste0(path,name_bir))
  HP = GenHpoly(i,3*i)
  tim=proc.time()
  st1=ban_volume(HP, steps_only = TRUE)
  tim=proc.time()-tim
  tim = as.numeric(as.character(tim[3]))
  time1=c(time1, tim)
  steps1=c(steps1, st1)
  
  tim=proc.time()
  st2=cg_volume(HP, steps_only = TRUE)
  tim=proc.time()-tim
  tim = as.numeric(as.character(tim[3]))
  time2=c(time2, tim)
  steps2=c(steps2, st2)
}

test.data <- data.frame(
  dimension= 2:100,
  BallAnnealing = steps1, 
  CoolingGaussian = steps2
)

save(test.data, file = "script_data.RData")

test.data %>%
  gather(methods,steps, BallAnnealing, CoolingGaussian) %>%
  ggplot(aes(x=dimension, y=steps, colour=methods)) +
  geom_line() + geom_point()