\name{lp.assign}
\alias{lp.assign}
\title{Integer Programming for the Assignment Problem}
\description{Interface to lp\_solve linear/integer programming
system specifically for solving assignment problems}
\usage{
lp.assign (cost.mat, direction = "min", presolve = 0, compute.sens = 0)
}
\arguments{
\item{cost.mat}{Matrix of costs: the ij-th element is the cost of 
assigning source i to destination j.}
\item{direction}{Character vector, length 1, containing either "min"
(the default) or "max"}
\item{presolve}{Numeric: presolve? Default 0 (no); any
non-zero value means "yes." Currently ignored.}
\item{compute.sens}{Numeric: compute sensitivity? Default 0 (no); any
non-zero value means "yes." In that case presolving is attempted.}
}
\details{
This is a particular integer programming problem. All the decision variables
are assumed to be integers; each row has the constraint that its entries 
must add up to 1 (so that there is one 1 and the remaining entries are 0) 
and each column has the same constraint. This is assumed to be a
minimization problem.
}
\value{
An \code{\link{lp}} object. See documentation for details. The constraints
are assumed (each row adds to 1, each column adds to 1, and no others) and
are not returned.
}
\author{Sam Buttrey, \email{buttrey@nps.edu}}
\seealso{\code{\link{lp}}, \code{\link{lp.transport}}}
\examples{
assign.costs <- matrix (c(2, 7, 7, 2, 7, 7, 3, 2, 7, 2, 8, 10, 1, 9, 8, 2), 4, 4)
\dontrun{
> assign.costs
     [,1] [,2] [,3] [,4]
[1,]    2    7    7    1
[2,]    7    7    2    9
[3,]    7    3    8    8
[4,]    2    2   10    2
}
lp.assign (assign.costs)
\dontrun{Success: the objective function is 8}
lp.assign (assign.costs)$solution
\dontrun{
     [,1] [,2] [,3] [,4]
[1,]    0    0    0    1
[2,]    0    0    1    0
[3,]    0    1    0    0
[4,]    1    0    0    0
}
}
\keyword{optimize}


