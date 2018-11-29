library(volesti)
#library(ggplot2)
#library(tidyr)

times1=c()
times2=c()
steps1=c()
steps2=c()
errors1=c()
errors2=c()
dimen=30
num_tests=1
path = system.file('extdata', package = 'volesti')
for (i in c(2,seq(from=5,to=40,by=5))){
  print(i)
  #name_bir = paste0('/birk',i,'.ine')
  #HP = fileToMatrix(paste0(path,name_bir))
  #HP = GenCube(i,'H')
  Z=GenZonotope(i,30*i)
  #ev=2^i
  
  st1=0
  err1=0
  tim1=0
  nb=0
  st2=0
  err2=0
  tim2=0
  #er11=c()
  for (j in 1:num_tests) {
    tim=system.time({ ps1 = ban_volume(Z)})
    tim = as.numeric(as.character(tim[3]))
    tim1=tim1+tim
    st1=st1+ps1[3]
    nb=ps1[2]
    #err1 = err1 + abs(ev-ps1[1])/ev
    #er11=c(er11,abs(ev-ps1[1])/ev)
    #print(paste0('vol = ',ps1[1]))
    
    #tim=system.time({ ps2 = cg_volume(HP)})
    #tim = as.numeric(as.character(tim[3]))
    #tim2=tim2+tim
    #st2=st2+ps2[3]
    #err2 = err2 + abs(ev-ps2[1])/ev
  }
  nb=nb/num_tests
  st1=st1/num_tests
  tim1=tim1/num_tests
  
  times1=c(times1,tim1)
  nballs1=c(nb,err1)
  steps1=c(steps1,st1)
  
  save(times1, file = "times_2_30_d_30d_ball.RData")
  save(steps1, file = "steps_2_30_d_30d_ball.RData")
  save(nballs1, file = "nballs_2_30_d_30d_ball.RData")
  
}

#test1.data <- data.frame(
#  dimension= 2:dimen,
#  CoolingBall = nballs1, 
#  CoolingGaussian = errors2
#


#save(test.data, file = "times_cubes_2_200.RData")
#save(test2.data, file = "steps_cubes_2_200.RData")

#test1.data %>%
#  gather(methods,steps, CoolingBall, CoolingGaussian) %>%
#  ggplot(aes(x=dimension, y=steps, colour=methods)) +
#  geom_line() + geom_point()
