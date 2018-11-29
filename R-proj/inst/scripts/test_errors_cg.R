library(volesti)
library(ggplot2)
library(tidyr)

times1=c()
times2=c()
steps1=c()
steps2=c()
errors1=c()
errors2=c()
dimen=200
num_tests=20
path = system.file('extdata', package = 'volesti')
for (i in 2:dimen) {
  print(i)
  #name_bir = paste0('/birk',i,'.ine')
  #HP = fileToMatrix(paste0(path,name_bir))
  HP = GenCube(i,'H')
  ev=2^i
  
  st1=0
  err1=0
  tim1=0
  st2=0
  err2=0
  tim2=0
  #er11=c()
  for (j in 1:num_tests) {
    #tim=system.time({ ps1 = ban_volume(HP)})
    #tim = as.numeric(as.character(tim[3]))
    #tim1=tim1+tim
    #st1=st1+ps1[3]
    #err1 = err1 + abs(ev-ps1[1])/ev
    #er11=c(er11,abs(ev-ps1[1])/ev)
    #print(paste0('vol = ',ps1[1]))
    
    tim=system.time({ ps2 = cg_volume(HP)})
    tim = as.numeric(as.character(tim[3]))
    tim2=tim2+tim
    st2=st2+ps2[3]
    err2 = err2 + abs(ev-ps2[1])/ev
  }
  err1=err1/num_tests
  st1=st1/num_tests
  tim1=tim1/num_tests
  
  err2=err2/num_tests
  st2=st2/num_tests
  tim2=tim2/num_tests
  
  times1=c(times1,tim1)
  errors1=c(errors1,err1)
  steps1=c(steps1,st1)
  
  times2=c(times2,tim2)
  errors2=c(errors2,err2)
  steps2=c(steps2,st2)
  
  save(times2, file = "times_cg_Hcubes_2_200.RData")
  save(steps2, file = "steps_cg_Hcubes_2_200.RData")
  save(errors2, file = "errors_cg_cubes_2_200.RData")
  
}

test1.data <- data.frame(
  dimension= 2:dimen,
  CoolingBall = errors1, 
  CoolingGaussian = errors2
)

test2.data <- data.frame(
  dimension= 2:dimen,
  CoolingBall = steps1, 
  CoolingGaussian = steps2
)

test3.data <- data.frame(
  dimension= 2:dimen,
  CoolingBall = times1, 
  CoolingGaussian = times2
)
#save(test.data, file = "times_cubes_2_200.RData")
#save(test2.data, file = "steps_cubes_2_200.RData")

test1.data %>%
  gather(methods,steps, CoolingBall, CoolingGaussian) %>%
  ggplot(aes(x=dimension, y=steps, colour=methods)) +
  geom_line() + geom_point()
