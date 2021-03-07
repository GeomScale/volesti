library(volesti)
P = gen_cube(2,'H')
p = sample_points(P,n=1000, distribution = list("density"="exponential",
                                               "bias"=c(0.5,0.5),"variance"=0.1))