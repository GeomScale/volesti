library(volesti)
library(ggplot2)
library(tidyr)

time1=c()
time2=c()
steps1=c()
steps2=c()
for (i in 2:40) {
  HP = GenCube(i,'H')
  tim=proc.time()
  ban_vec=ban_volume(HP)
  tim=proc.time()-tim
  tim = as.numeric(as.character(tim[3]))
  time1=c(time1, tim)
  steps1=c(steps1, ban_vec[3])
  
  tim=proc.time()
  cg_vec=cg_volume(HP)
  tim=proc.time()-tim
  tim = as.numeric(as.character(tim[3]))
  time2=c(time2, tim)
  steps2=c(steps2, cg_vec[3])
}

test.data <- data.frame(
  dimension= 2:40,
  ban_steps = steps1, 
  cg_steps = steps2
)

test.data %>%
  gather(methods,steps, ban_steps, cg_steps) %>%
  ggplot(aes(x=dimension, y=steps, colour=methods)) +
  geom_line() + geom_point()