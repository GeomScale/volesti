
library(volesti)
rets = opti_hmc(50,50,70,15000)

save(rets, file = "sdp_hmc_hnr3.RData")



boundcalls_vpoly.data <- data.frame(
  value = c(rets[1,], rets[3,]),
  random_walk = c(rep("HMC",70), rep("RDHR",70)),
  iteration = rep(1:70,2)
)

ggplot(boundcalls_vpoly.data , aes(x=iteration, y=value)) +
  geom_line(aes(color=random_walk),size=1) + geom_point() +labs(y="Obj. function values")  +
  scale_y_continuous(breaks = c(-1.3,-1.209, -1, -0.5, 0)) + scale_x_continuous(breaks = seq(from=0,to=70,by=10)) +
  geom_hline(yintercept=-1.209, linetype="dashed", color = "black", size=0.4) + labs(color="random walk")+
  theme(legend.position="top",text = element_text(size=30), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))

