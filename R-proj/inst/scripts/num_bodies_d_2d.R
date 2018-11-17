library(volesti)
#library(ggplot2)
#library(tidyr)

bodies1=c()
#bodies2=c()
continue_ball = TRUE
time_limit = 2700
i=80
while(TRUE) {
  print(paste0('i = ',i))
  if (!continue_ball) {
    break
  }
  num_b=0
  #num_b2 = 0
  for (j in 1:1){
    Z = GenZonotope(i,2*i)
    x=system.time({ b1 = ban_volume(Z,steps_only = TRUE)})
    if (as.numeric(x[3])>time_limit) {
      continue_ball = TRUE
    }
    num_b = num_b + b1
  }
    
  num_b = num_b / 1
  bodies1=c(bodies1, num_b)
  save(bodies1, file = "bodies80_d_2d.RData")
  if(i==2){
    i=5
  }else{
    i = i + 5
  }
}

