library(volesti)
library(ggplot2)
library(tidyr)

errors1 = c()
steps1 = c()
errors2 = c()
steps2 = c()
for (i in 2:100){
  
  #print(paste0('i = ', i))
  #H=GenCube(i,'H')
  H=GenHpoly(i,4*i)
  ex_vol = 2^i
  vol_it1 = 0
  steps_it1 = 0
  vol_it2 = 0
  steps_it2 = 0
  for (j in 1:5){
    st=ban_volume(H,len_subwin = 30, len_tuple = 150+i)
    #st2=cg_volume(H)
    st2=ban_volume(H,const_win = FALSE)
    vol_it1 = vol_it1 + st[1]
    steps_it1 = steps_it1 + st[3]
    vol_it2 = vol_it2 + st2[1]
    steps_it2 = steps_it2 + st2[3]
  }
  vol_it1 = vol_it1 / 5
  steps_it1 = steps_it1 / 5
  vol_it2 = vol_it2 / 5
  steps_it2 = steps_it2 / 5
  
  err = abs(ex_vol - vol_it1)/ex_vol
  errors1 = c(errors1, err)
  steps1 = c(steps1, steps_it1)
    
  err = abs(ex_vol - vol_it2)/ex_vol
  errors2 = c(errors2, err)
  steps2 = c(steps2, steps_it2)
  
}

test1.data <- data.frame(
  dimension= 2:(i),
  new_window = errors1, 
  Cousins_window = errors2
)

test1.data %>%
  gather(methods,errors,new_window, Cousins_window) %>%
  ggplot(aes(x=dimension, y=errors, colour=methods)) +
  geom_line() + geom_point()

test2.data <- data.frame(
  dimension= 2:(i),
  new_window = steps1, 
  Cousins_window = steps2
)

#save(test2.data, file = "d_15d_data.RData")

test2.data %>%
  gather(methods,steps,new_window, Cousins_window) %>%
  ggplot(aes(x=dimension, y=steps, colour=methods)) +
  geom_line() + geom_point()
