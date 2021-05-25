\name{print.lp}
\alias{print.lp}
\title{Print an lp object}
\description{Print method for lp objects}
\usage{
\method{print}{lp} (x, \ldots)
}
\arguments{
\item{x}{List with items named \code{objval} and \code{status}.
Normally this will have been called by \code{\link{lp}},
\code{\link{lp.assign}}, or \code{\link{lp.transport}}.}
\item{...}{Other arguments, all currently ignored}
}
\details{
This function prints the objective function value, together with the
word "Success" if the operation is successful, or an indication of the
error if not. If multiple solutions have been produced (because this was
an all-binary problem and lp was called with num.bin.solns > 1) the number
of solutions is also displayed.}
\value{
None
}
\author{Sam Buttrey, \email{buttrey@nps.edu}}
\seealso{\code{\link{lp}}, \code{\link{lp.assign}}, \code{\link{lp.transport}}}
\keyword{optimize}
