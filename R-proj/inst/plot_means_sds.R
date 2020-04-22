library(ggplot2)
library(matrixStats)

load("~/volume_approximation/R-proj/inst/ps.RData")

frts_rows = c(1,4,7,10,13,16,19,22,25,28,31)
first_bill = ps[frts_rows,]
second_bill = ps[frts_rows+1,]
third_bill = ps[frts_rows+2,]

frts_rows = c(34,37,40,43,46,49,52,55,58,61,64)
first_rdhr = ps[frts_rows,]
second_rdhr = ps[frts_rows+1,]
third_rdhr = ps[frts_rows+2,]

frts_rows = c(67,70,73,76,79,82,85,88,91,94,97)
first_cdhr = ps[frts_rows,]
second_cdhr = ps[frts_rows+1,]
third_cdhr = ps[frts_rows+2,]

boundcalls_vpoly.data <- data.frame(
  Mean = c(rowMeans(first_bill), rowMeans(first_rdhr), rowMeans(first_cdhr),rowMeans(second_bill), rowMeans(second_rdhr), rowMeans(second_cdhr), rowMeans(third_bill), rowMeans(third_rdhr), rowMeans(third_cdhr)),
  random_walk = c(rep("BW",11), rep("RDHR",11),rep("CDHR",11), rep("BW",11), rep("RDHR",11),rep("CDHR",11), rep("BW",11), rep("RDHR",11),rep("CDHR",11)),
  ratio = c(rep("0.0993",11), rep("0.0993",11), rep("0.0993",11), rep("0.588",11), rep("0.588",11), rep("0.588",11), rep("0.716",11), rep("0.716",11), rep("0.716",11)),
  walk_length = rep(c(1,seq(from=5,to=50,by=5)),3)
)

ggplot(boundcalls_vpoly.data , aes(x=walk_length, y=Mean)) +
       geom_line(aes(linetype =ratio, color=random_walk),size=1) + geom_point() +labs(x ="walk length", y=expression(paste("unbiased estimator ", hat(E),"[f]")))  +
       scale_y_continuous(breaks = c(0,0.0993,0.588,0.716)) + geom_hline(yintercept=0.0993, linetype="dashed", color = "black", size=0.4) +
       geom_hline(yintercept=0.5880, linetype="dashed", color = "black", size=0.4) + geom_hline(yintercept=0.7160, linetype="dashed", color = "black", size=0.4) + guides(fill=guide_legend(title="New Legend Title"))+
       labs(linetype=expression(paste(math(E),"[f]")), color="random walk")+
       theme(legend.position="top",text = element_text(size=30), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))



boundcalls_vpoly.data <- data.frame(
  Standard_deviation = c(rowSds(first_bill), rowSds(first_rdhr), rowSds(first_cdhr),rowSds(second_bill), rowSds(second_rdhr), rowSds(second_cdhr), rowSds(third_bill), rowSds(third_rdhr), rowSds(third_cdhr)),
  random_walk = c(rep("BW",11), rep("RDHR",11),rep("CDHR",11), rep("BW",11), rep("RDHR",11),rep("CDHR",11), rep("BW",11), rep("RDHR",11),rep("CDHR",11)),
  ratio = c(rep("0.0993",11), rep("0.0993",11), rep("0.0993",11), rep("0.588",11), rep("0.588",11), rep("0.588",11), rep("0.716",11), rep("0.716",11), rep("0.716",11)),
  walk_length = rep(c(1,seq(from=5,to=50,by=5)),3)
)

 ggplot(boundcalls_vpoly.data , aes(x=walk_length, y=Standard_deviation)) +
       geom_line(aes(linetype =ratio, color=random_walk),size=1) + geom_point() +labs(x ="walk length",y=expression(paste("st.d. of unbiased estimator ",hat(E),"[f]")))  +
       labs(linetype=expression(paste(math(E),"[f]")), color="random walk")+
       theme(legend.position="top",text = element_text(size=30), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))
 
 
 #-----------------------only two ratios ----------------------#
 
 cbp1 <- c("red", "blue", "gold")
 lt = c("solid", "dotdash")
 
 boundcalls_vpoly.data <- data.frame(
   Mean = c(rowMeans(first_bill), rowMeans(first_rdhr), rowMeans(first_cdhr),rowMeans(second_bill), rowMeans(second_rdhr), rowMeans(second_cdhr)),
   random_walk = c(rep("W-Billiard",11), rep("W-HNR",11),rep("W-CHNR",11), rep("W-Billiard",11), rep("W-HNR",11),rep("W-CHNR",11)),
   ratio = c(rep("0.0993",11), rep("0.0993",11), rep("0.0993",11), rep("0.588",11), rep("0.588",11), rep("0.588",11)),
   walk_length = rep(c(1,seq(from=5,to=50,by=5)),6)
 )
 
 ggplot(boundcalls_vpoly.data , aes(x=walk_length, y=Mean)) +
   geom_line(aes(linetype =ratio, color=random_walk),size=1) + geom_point() +labs(x ="walk length", y=expression(paste("unbiased estimator ", hat(E),"[f]")))  +
   scale_y_continuous(breaks = c(0,0.0993,0.588)) + geom_hline(yintercept=0.0993, linetype="dashed", color = "black", size=0.4) +
   geom_hline(yintercept=0.5880, linetype="dashed", color = "black", size=0.4) +
   labs(linetype=expression(paste(math(E),"[f]")), color="random walk")+ scale_color_manual(values = cbp1) +scale_linetype_manual(values = lt) +
   theme(legend.position="top",text = element_text(size=30), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))
 
 
 
 
 boundcalls_vpoly.data <- data.frame(
   Standard_deviation = c(rowSds(first_bill), rowSds(first_rdhr), rowSds(first_cdhr),rowSds(second_bill), rowSds(second_rdhr), rowSds(second_cdhr)),
   random_walk = c(rep("W-Billiard",11), rep("W-HNR",11),rep("W-CHNR",11), rep("W-Billiard",11), rep("W-HNR",11),rep("W-CHNR",11)),
   ratio = c(rep("0.0993",11), rep("0.0993",11), rep("0.0993",11), rep("0.588",11), rep("0.588",11), rep("0.588",11)),
   walk_length = rep(c(1,seq(from=5,to=50,by=5)),6)
 )
 
 ggplot(boundcalls_vpoly.data , aes(x=walk_length, y=Standard_deviation)) +
   geom_line(aes(linetype =ratio, color=random_walk),size=1) + geom_point() +labs(x ="walk length",y=expression(paste("st.d. of unbiased estimator ",hat(E),"[f]")))  +
   labs(linetype=expression(paste(math(E),"[f]")), color="random walk")+scale_color_manual(values = cbp1) +scale_linetype_manual(values = lt) +
   theme(legend.position="top",text = element_text(size=30), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))
 
 
 