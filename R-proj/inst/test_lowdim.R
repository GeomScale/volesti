library(volesti)

P = GenSimplex(3,'H')

A = P$A[1:3,]
b = P$b

Aeq = matrix(c(1, 1, 1), nrow=1, ncol=3)
beq = c(1)

p = SamplePoints(A=A, b=b, Aeq = Aeq, beq = beq, n=800, boundary = TRUE , random_walk = "RDHR", walk_length = 10)