library(volesti)
library(ggplot2)
library(tidyr)

bodies1=c()
bodies2=c()
for (i in 2:5) {
  Z = GenZonotope(i,2*i)
  num_b=ban_volume(Z,steps_only = TRUE)
  bodies1=c(bodies1, num_b)
  

  num_b=cg_volume(Z,steps_only = TRUE)
  bodies2=c(bodies2, num_b)
}

test2.data <- data.frame(
  dimension= 2:5,
  ban_bodies = bodies1, 
  cg_bodies = bodies2
)

test2.data %>%
  gather(methods,bodies,ban_bodies, cg_bodies) %>%
  ggplot(aes(x=dimension, y=bodies, colour=methods)) +
  geom_line() + geom_point()