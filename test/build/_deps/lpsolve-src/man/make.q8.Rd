\name{make.q8}
\alias{make.q8}
\alias{8-queens problem}
\title{Generate sparse constraint matrix for 8-queens problem}
\description{Generate sparse constraint matrix for 8-queens problem}
\usage{ make.q8 () }
\arguments{None.}
\details{
Sparse constraints come in a three-column matrix or data frame. Each row
gives the row number, column number, and value of a particular non-zero
entry in the constraint matrix. This function produces the sparse constraint
matrix for the 8-queens problem (in which the object is to place eight queens
on a chessboard with no two sharing a row, column or diagonal). The resulting
sparse representation is 252 x 3, compared to 42 x 64 for the usual 
representation.}
\value{
A 252 x 3 numeric matrix. See \link{lp} for the complete example.
}
\author{Sam Buttrey, \email{buttrey@nps.edu}}
\seealso{\code{\link{lp}}}
\keyword{optimize}
