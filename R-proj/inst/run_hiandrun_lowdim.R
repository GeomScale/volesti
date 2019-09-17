 P = GenSimplex(3,'H')
 ineq.constr <- list(
       constr = P$A[1:3,],
       dir = rep('<=', 3),
       rhs = c(rep(0, 3)))
eq.constr <- list(
       constr = matrix(c(1, 1, 1), nrow=1, ncol=3),
       dir = '=',
       rhs = 1)

basis <- hitandrun::solution.basis(eq.constr)
transform <- hitandrun::createTransform(basis)
constr <- hitandrun::transformConstraints(transform, ineq.constr)
x0 <- hitandrun::createSeedPoint(constr, homogeneous=TRUE)
x <- hitandrun::har(x0, constr, 500, transform=transform, homogeneous=TRUE)$samples