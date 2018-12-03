library(volesti)
library(geometry)
#library(ggplot2)
#library(tidyr)

times1=c()
times2=c()
steps1=c()
steps2=c()
errors1=c()
errors2=c()
nballs1=c()
vols1=c()
dimen=100
num_tests=1
path = system.file('extdata', package = 'volesti')
for (i in c(seq(from=50,to=100,by=5))) {
  print(i)
  #name_bir = paste0('/birk',i,'.ine')
  #HP = fileToMatrix(paste0(path,name_bir))
  #HP = GenCube(i,'H')
  
  #ev=2^i/prod(1:i)
  
  st1=0
  #err1=0
  tim1=0
  vol1=0
  nb=0
  #st2=0
  #err2=0
  #tim2=0
  #er11=c()
  for (j in 1:num_tests) {
    P=GenVpoly(i,10*i)
    
    tim=system.time({ ps1 = ban_volume(P,rounding = TRUE)})
    tim = as.numeric(as.character(tim[3]))
    tim1=tim1+tim
    st1=st1+ps1[3]
    nb=ps1[2]+nb
    vol1=vol1+ps1[1]
    #err1 = err1 + abs(ev-ps1[1])/ev
    #er11=c(er11,abs(ev-ps1[1])/ev)
    #print(paste0('vol = ',ps1[1]))
    
    #tim=system.time({ ps2 = cg_volume(HP)})
    #tim = as.numeric(as.character(tim[3]))
    #tim2=tim2+as.numeric(as.character(tim_ex[3]))
    #st2=st2+ps2[3]
    #err2 = err2 + abs(ev-ps2[1])/ev
  }
  vol1=vol1/num_tests
  nb=nb/num_tests
  st1=st1/num_tests
  tim1=tim1/num_tests
  #tim2=tim2/num_tests
  #err1 = err1/num_tests
  
  times1=c(times1,tim1)
  vols1 =c(vols1,vol1)
  nballs1=c(nballs1,nb)
  steps1=c(steps1,st1)
  #times2=c(times2,tim2)
  #errors1=c(errors1,err1)
  
  save(times1, file = "1xtimes1_50_100_randV_d_10d.RData")
  #save(times2, file = "times2_2_14_randV_d_6d.RData")
  #save(errors1, file = "errors_2_14_randV_d_6d.RData")
  save(vols1, file = "1xvols_50_100_randV_d_10d.RData")
  save(steps1, file = "1xsteps_50_100_randV_d_10d.RData")
  save(nballs1, file = "1xnballs_50_100_randV_d_10d.RData")
  
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
