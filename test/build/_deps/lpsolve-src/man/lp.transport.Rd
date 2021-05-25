\name{lp.transport}
\alias{lp.transport}
\title{Integer Programming for the Transportation Problem}
\description{Interface to lp\_solve linear/integer programming
system specifically for solving transportation problems}
\usage{
lp.transport (cost.mat, direction="min", row.signs, row.rhs, col.signs,
              col.rhs, presolve=0, compute.sens=0, integers = 1:(nc*nr) )
}
\arguments{
\item{cost.mat}{Matrix of costs; ij-th element is the cost of transporting 
one item from source i to destination j.}
\item{direction}{Character, length 1: "min" or "max"}
\item{row.signs}{Vector of character strings giving the direction of the 
row constraints: each value should be one of "<," "<=," "=," "==," ">," 
or ">=." (In each pair the two values are identical.)}
\item{row.rhs}{Vector of numeric values for the right-hand sides of the 
row constraints.}
\item{col.signs}{Vector of character strings giving the direction of the 
column constraints: each value should be one of "<," "<=," "=," "==," ">," 
or ">=."}
\item{col.rhs}{Vector of numeric values for the right-hand sides of the 
column constraints.}
\item{presolve}{Numeric: presolve? Default 0 (no); any
non-zero value means "yes." Currently ignored.}
\item{compute.sens}{Numeric: compute sensitivity? Default 0 (no); any
non-zero value means "yes."}
\item{integers}{Vector of integers whose ith element gives the index 
of the ith integer variable. Its length will be the number of integer 
variables. Default: all variables are integer. Set to NULL to have no
variables be integer.}
}
\details{
This is a particular integer programming problem. All the decision variables
are assumed to be integers, and there is one constraint per row and one per
column (and no others). This is assumed to be a minimization problem.
}
\value{
An \code{\link{lp}} object. Constraints are implicit and not returned.
See documentation for details.
}
\references{Example problem from Bronson (1981), \emph{Operations Research},
Scahum's Outline Series, McGraw-Hill.}
\author{Sam Buttrey, \email{buttrey@nps.edu}}
\seealso{\code{\link{lp.assign}}, \code{\link{lp.transport}}}

\examples{
#
# Transportation problem, Bronson, problem 9.1, p. 86
#
# Set up cost matrix
#
costs <- matrix (10000, 8, 5); costs[4,1] <- costs[-4,5] <- 0
costs[1,2] <- costs[2,3] <- costs[3,4] <- 7; costs[1,3] <- costs[2,4] <- 7.7
costs[5,1] <- costs[7,3] <- 8; costs[1,4] <- 8.4; costs[6,2] <- 9
costs[8,4] <- 10; costs[4,2:4] <- c(.7, 1.4, 2.1)
#
# Set up constraint signs and right-hand sides.
#
row.signs <- rep ("<", 8)
row.rhs <- c(200, 300, 350, 200, 100, 50, 100, 150)
col.signs <- rep (">", 5)
col.rhs <- c(250, 100, 400, 500, 200)
#
# Run
#
lp.transport (costs, "min", row.signs, row.rhs, col.signs, col.rhs)
\dontrun{Success: the objective function is 7790}
lp.transport (costs, "min", row.signs, row.rhs, col.signs, col.rhs)$solution
\dontrun{
     [,1] [,2] [,3] [,4] [,5]
[1,]    0  100    0  100    0
[2,]    0    0  300    0    0
[3,]    0    0    0  350    0
[4,]  200    0    0    0    0
[5,]   50    0    0    0   50
[6,]    0    0    0    0   50
[7,]    0    0  100    0    0
[8,]    0    0    0   50  100
}
}
\keyword{optimize}


