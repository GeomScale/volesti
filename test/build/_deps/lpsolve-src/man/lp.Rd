\name{lp}
\alias{lp}
\title{Linear and Integer Programming}
\description{Interface to lp\_solve linear/integer programming system}
\usage{
lp (direction = "min", objective.in, const.mat, const.dir, const.rhs,
	transpose.constraints = TRUE, int.vec, presolve=0, compute.sens=0,
        binary.vec, all.int=FALSE, all.bin=FALSE, scale = 196, dense.const, 
        num.bin.solns=1, use.rw=FALSE)
}
\arguments{
\item{direction}{Character string giving direction of optimization:
"min" (default) or "max."}
\item{objective.in}{Numeric vector of coefficients of objective function}
\item{const.mat}{Matrix of numeric constraint coefficients, one row 
per constraint, one column per variable (unless transpose.constraints = 
FALSE; see below).}
\item{const.dir}{Vector of character strings giving the direction of 
the constraint: each value should be one of "<," "<=," "=," "==," ">," or ">=".
(In each pair the two values are identical.)}
\item{const.rhs}{Vector of numeric values for the right-hand sides of 
the constraints.}
\item{transpose.constraints}{By default each constraint occupies a row 
of const.mat, and that matrix needs to be transposed before being passed 
to the optimizing code.  For very large constraint matrices it may be wiser 
to construct the constraints in a matrix column-by-column. In that case set 
transpose.constraints to FALSE.}
\item{int.vec}{Numeric vector giving the indices of variables that are 
required to be integer. The length of this vector will therefore be the 
number of integer variables.}
\item{presolve}{Numeric: presolve? Default 0 (no); any
non-zero value means "yes." Currently ignored.}
\item{compute.sens}{Numeric: compute sensitivity? Default 0 (no); any
non-zero value means "yes."}
\item{binary.vec}{Numeric vector like int.vec giving the indices of variables 
that are required to be binary.}
\item{all.int}{Logical: should all variables be integer? Default: FALSE.}
\item{all.bin}{Logical: should all variables be binary? Default: FALSE.}
\item{scale}{Integer: value for lpSolve scaling. Details can be found in
the lpSolve documentation. Set to 0 for no scaling. Default: 196}
\item{dense.const}{Three column dense constraint array. This is ignored if
const.mat is supplied. Otherwise the columns are constraint number, column
number, and value; there should be one row for each non-zero entry in the
constraint matrix.}
\item{num.bin.solns}{Integer: if all.bin=TRUE, the user can request up to
num.bin.solns optimal solutions to be returned.}
\item{use.rw}{Logical: if TRUE and num.bin.solns > 1, write the lp out to a
file and read it back in for each solution after the first. This is just to
defeat a bug somewhere. Although the default is FALSE, we recommend you set
this to TRUE if you need num.bin.solns > 1, until the bug is found.}
}
\details{
This function calls the lp\_solve 5.5 solver. That system has many options not
supported here. The current version is maintained at
\url{http://lpsolve.sourceforge.net/5.5}

Note that every variable is assumed to be >= 0!
}
\value{
An lp object. See \code{\link{lp.object}} for details.
}
\author{Sam Buttrey, \email{buttrey@nps.edu}}
\seealso{\code{\link{lp.assign}}, \code{\link{lp.transport}}}
\examples{
#
# Set up problem: maximize
#   x1 + 9 x2 +   x3 subject to
#   x1 + 2 x2 + 3 x3  <= 9
# 3 x1 + 2 x2 + 2 x3 <= 15
#
f.obj <- c(1, 9, 1)
f.con <- matrix (c(1, 2, 3, 3, 2, 2), nrow=2, byrow=TRUE)
f.dir <- c("<=", "<=")
f.rhs <- c(9, 15)
#
# Now run.
#
lp ("max", f.obj, f.con, f.dir, f.rhs)
\dontrun{Success: the objective function is 40.5}
lp ("max", f.obj, f.con, f.dir, f.rhs)$solution
\dontrun{[1] 0.0 4.5 0.0}
#
# The same problem using the dense constraint approach:
#
f.con.d <- matrix (c(rep (1:2,each=3), rep (1:3, 2), t(f.con)), ncol=3)
lp ("max", f.obj, , f.dir, f.rhs, dense.const=f.con.d)
\dontrun{Success: the objective function is 40.5}
#
# Get sensitivities
#
lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$sens.coef.from
\dontrun{[1] -1e+30  2e+00 -1e+30}
lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$sens.coef.to  
\dontrun{[1] 4.50e+00 1.00e+30 1.35e+01}
#
# Right now the dual values for the constraints and the variables are
# combined, constraints coming first. So in this example...
#
lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$duals     
\dontrun{[1]   4.5   0.0  -3.5   0.0 -10.5}
#
# ...the duals of the constraints are 4.5 and 0, and of the variables,
# -3.5, 0.0, -10.5. Here are the lower and upper limits on these:
#
lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$duals.from
\dontrun{[1]  0e+00 -1e+30 -1e+30 -1e+30 -6e+00}
lp ("max", f.obj, f.con, f.dir, f.rhs, compute.sens=TRUE)$duals.to  
\dontrun{[1] 1.5e+01 1.0e+30 3.0e+00 1.0e+30 3.0e+00}
#
# Run again, this time requiring that all three variables be integer
#
lp ("max", f.obj, f.con, f.dir, f.rhs, int.vec=1:3)
\dontrun{Success: the objective function is 37}
lp ("max", f.obj, f.con, f.dir, f.rhs, int.vec=1:3)$solution
\dontrun{[1] 1 4 0}
#
# You can get sensitivities in the integer case, but they're harder to
# interpret.
#
lp ("max", f.obj, f.con, f.dir, f.rhs, int.vec=1:3, compute.sens=TRUE)$duals
\dontrun{[1] 1 0 0 7 0}
#
# Here's an example in which we want more than one solution to a problem
# in which all variables are binary: the 8-queens problem, 
# with dense constraints.
#
chess.obj <- rep (1, 64)
q8 <- make.q8 ()
chess.dir <- rep (c("=", "<"), c(16, 26))
chess.rhs <- rep (1, 42)
lp ('max', chess.obj, , chess.dir, chess.rhs, dense.const = q8, 
    all.bin=TRUE, num.bin.solns=3)
}
\keyword{optimize}
