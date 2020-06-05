\name{Spectrahedron}
\alias{Spectrahedron}
\title{An \code{R} class to represent spectrahedra.}

\description{
A spectrahedron is a convex body defined by a linear matrix inequality of the form \eqn{A_0 + x_1 A_1 + ... + x_n A_n \preceq 0}.
The matrices \eqn{A_i} are symmetric \eqn{m \times m} real matrices and \eqn{\preceq 0} denoted negative semidefiniteness.
}
\section{Fields}{
\itemize{
\item{\code{matrices} }{A list with the matrices \eqn{A_0, A_1, ..., A_n}}
}}
