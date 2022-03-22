# A comparative study of uniform high dimensional samplers

## Sampling in a polytope.
```
walk="CDHR"
b=c(10,10,10,10,10)
A = matrix(c(1,0,-0.25,-1,2.5,1,0.4,-1,-0.9,0.5), nrow=5, ncol=2, byrow = TRUE)
P = Hpolytope(A = A, b = b)
points = sample_points(P, 1000)
plot(ggplot(data.frame( x=points[1,], y=points[2,] )) +
geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-15,15),
ylim = c(-15,15)) + ggtitle(sprintf("Sampling a random pentagon with walk %s", walk)))
```
![Sampling a polygon.](/images/pentagon.png "Sampling a polytope.")

## Compair random walks on the 100-dimensional hypercube.
After sampling 1000 points for every walk we project on the first two coordinates and essentially comment on their mixing time. 
```
for (step in c(1,20,50,100,150)){
  for (walk in c("CDHR", "RDHR", "BaW")){
    P <- gen_cube(100, 'H')
    points1 <- sample_points(P, 1000, random_walk = list("walk" = walk, "walk_length" = step))
    jpeg(file=paste(step,walk,'.jpg'));
    g<-plot(ggplot(data.frame( x=points1[1,], y=points1[2,] )) +
geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-1,1),
ylim = c(-1,1)) + ggtitle(sprintf("walk length=%s", step)))
dev.off()
  }
}
```
For walk length 1 we can see that neither one of the random walk sampler has reached the mixing time.
![Sampling a polygon.](/images/1 CDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/1 RDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/1 BaW.jpg "Sampling a polytope.")

For walk length 20 we can see that the CDHR random sampler seems to be closer to the mixing time. While the other walks are far from it.

![Sampling a polygon.](/images/20 CDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/20 RDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/20 BaW.jpg "Sampling a polytope.")

For walk length 50 we can see that the CDHR random sampler seems to have reached the mixing time. While the other walks still have not.

![Sampling a polygon.](/images/50 CDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/50 RDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/50 BaW.jpg "Sampling a polytope.")

For walk length 100 we can see that the CDHR and RDHR random sampler has definately reached the mixing time. While the BaW has not yet explored the whole space.

![Sampling a polygon.](/images/100 CDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/100 RDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/100 BaW.jpg "Sampling a polytope.")

For walk length 150 we can see that all of the samplers seem to produce uniformly distributed points.
![Sampling a polygon.](/images/150 CDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/150 RDHR.jpg "Sampling a polytope.")
![Sampling a polygon.](/images/150 BaW.jpg "Sampling a polytope.")

So the mixing time we can postulate that the mixing time for CDHR is between 20 and 50, mixing time of RDHR is between 50 and 100 and the mixing time of BaW is between 100 and 150.




## Statistical test for convergence of a random walk sampler.
We will use the Gelman-Rubin test for MCMC convergence. This test when given multiple MCMC chains compares in chain variability to the between chain variability.  As all of the coordinates need to be independent and converge to the same distribution we can apply it between them.

```
library(coda)
for (step in c(1,20,50,100,150)){
  for (walk in c("CDHR", "RDHR", "BaW")){
    P <- gen_cube(100, 'H')
    points1 <- sample_points(P, 1000, random_walk = list("walk" = walk, "walk_length" = step))
    jpeg(file=paste(step,walk,'.jpg'));
    g<-plot(ggplot(data.frame( x=points1[1,], y=points1[2,] )) +
geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-1,1),
ylim = c(-1,1)) + ggtitle(sprintf("walk length=%s", step)))
dev.off()
print(sprintf("walk_type=%s walk_length=%s",walk,step))
a=as.mcmc.list(lapply(as.data.frame(t(points1)), mcmc));
print(gelman.diag(a,autoburnin = FALSE,multivariate = TRUE))
jpeg(file=paste(step,walk,'gel','.jpg'));
gelman.plot(a,autoburnin = FALSE)
dev.off()
  }

```
The results for the above walks are:

| Walk | Walk length | R |
| :---: | :---: | :---: |
| CDHR | 1 |1.11  |
| RDHR | 1 | 1.55 |
| BaW | 1 | 1.99 |
| CDHR | 20 |  1 |
| RDHR | 20 |  1.06 |
| BaW | 20| 1.11 |
| CDHR | 50 |  1 |
| RDHR | 50 |  1.02 |
| BaW | 50 | 1.05 |
| CDHR | 100 |  1 |
| RDHR | 100 |  1.01 |
| BaW | 100 | 1.02 |
| CDHR | 150 |  1 |
| RDHR | 150 |  1.01 |
| BaW | 150 | 1.02 |
We know that the the walk is close to convergence if the value of the Gelman-Rubin statistic is close to one. From this statistic we can see that we are close to convergence for  walk length greater than equal to 20 for CDHR and RDHR. And BaW is close to convergence for walk length greater than equal to 50.
