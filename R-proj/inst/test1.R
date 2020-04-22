library(ggplot2)
library(volesti)

ps = spectra_plot(2, 7, 5000, 800, 5)
 points1 = ps[1:2,1:10000]
points3 = ps[1:2,10801:11600]
points2 = ps[1:2,10001:10800]
 
 cbp1 <- c("grey3", "red4")
 ggplot(data.frame(x = c(points3[1,],points1[1,]), body = c(rep("P",800),
         rep("boundary",10000)), y = c(points3[2,],points1[2,])) , aes(x=x, y=y)) +
         scale_color_manual(values =cbp1) +
        geom_point(shape=20, aes(color=body)) +labs(x =" ", y = " ")
 
 uniform_points = sample_points('sdp_prob_2_6.txt', N=2000)
 
 ggplot(data.frame(x = uniform_points[1,], y = uniform_points[2,]),aes(x=x, y=y)) +
       geom_point(shape=20,color="red") +labs(x =" ", y = " ")+xlim(-1.6, 1.8)+ylim(-2,2.3)
 
 boltz_points = sample_points('sdp_prob_2_6.txt', N=2000, distribution = 'boltzmann', Temperature = 2)
 
 ggplot(data.frame(x = boltz_points[1,], y = boltz_points[2,]),aes(x=x, y=y)) +
   geom_point(shape=20,color="red") +labs(x =" ", y = " ")+xlim(-1.6, 1.8)+ylim(-2,2.3)
 
 plot_ly(x = ~uniform_points[1,], y = ~uniform_points[2,], z = ~uniform_points[3,],marker = list(color = 'red', colorscale = c('red'),size=1, showscale = TRUE))
 plot_ly(x = ~boltz_points[1,], y = ~boltz_points[2,], z = ~boltz_points[3,],marker = list(color = 'red', colorscale = c('red'),size=1, showscale = TRUE))
 