\name{lp.object}
\alias{lp.object}
\title{LP (linear programming) object}
\description{Structure of lp object}
\value{
An lp.object is a list containing the following elements:
\item{direction}{Optimization direction, as entered}
\item{x.count}{Number of variables in objective function}
\item{objective}{Vector of objective function coefficients, as entered}
\item{const.count}{Number of constraints entered}
\item{constraints}{Constraint matrix, as entered (not returned
by \code{\link{lp.assign}} or \code{\link{lp.transport}})}
\item{int.count}{Number of integer variables}
\item{int.vec}{Vector of integer variables' indices, as entered}
\item{objval}{{Value of objective function at optimum}}
\item{solution}{Vector of optimal coefficients}
\item{num.bin.solns}{Numeric indicator of number of solutions returned}
\item{status}{Numeric indicator: 0 = success, 2 = no feasible solution}
}
\author{Sam Buttrey, \email{buttrey@nps.edu}}
\seealso{\code{\link{lp}}, \code{\link{lp.assign}}, \code{\link{lp.transport}}}
\keyword{optimize}
