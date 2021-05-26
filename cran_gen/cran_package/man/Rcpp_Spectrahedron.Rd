\docType{class}
\name{Rcpp_Spectrahedron}
\alias{Rcpp_Spectrahedron-class}
\alias{[,Rcpp_Spectrahedron-method}
\alias{[,Rcpp_Spectrahedron,ANY,ANY,ANY-method}
\alias{$<-,Rcpp_Spectrahedron-method}
\alias{$,Rcpp_Spectrahedron-method}
\alias{filepaths<-,Rcpp_Spectrahedron-method}
\title{
An \code{Rcpp} class to represent spectrahedra, exposed to \code{R} via modules.
}
\description{
A spectrahedron is a convex body defined by a linear matrix inequality of the form \eqn{A_0 + x_1 A_1 + ... + x_n A_n \preceq 0}.
The matrices \eqn{A_i} are symmetric \eqn{m \times m} real matrices and \eqn{\preceq 0} denoted negative semidefiniteness.
}
\details{
\describe{
\item{\code{matrices} }{A list with the matrices \eqn{A_0, A_1, ..., A_n}}
 }
}

