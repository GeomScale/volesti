\docType{class}
\name{Rcpp_PolytopeIntersectEllipsoid}
\alias{Rcpp_PolytopeIntersectEllipsoid-class}
\alias{[,Rcpp_PolytopeIntersectEllipsoid-method}
\alias{[,Rcpp_PolytopeIntersectEllipsoid,ANY,ANY,ANY-method}
\alias{$<-,Rcpp_PolytopeIntersectEllipsoid-method}
\alias{$,Rcpp_PolytopeIntersectEllipsoid-method}
\alias{filepaths<-,Rcpp_PolytopeIntersectEllipsoid-method}
\title{
An \code{Rcpp} class to represent an intersection between an H-polytope and an ellipsoid, exposed to \code{R} via modules.
}
\description{
An intersection between of an H-polytope, defined by a set of linear inequalities, or a matrix A and an \eqn{m}-dimensional vector b s.t., \eqn{x, Ax\leq b}, and an ellipsoid defined by a positive definite matrix E s.t., \eqn{x, x^TEx \leq 1}. 
}
\details{
\describe{
\item{\code{A} }{ \eqn{m \times d} matrix A}

\item{\code{b} }{ \eqn{m}-dimensional vector b}

\item{\code{E} }{ \eqn{d \times d} positive definite matrix E}

\item{\code{type} }{ An integer that declares the representation of the convex body. For the intersection between an H-polytope and an ellipsoid the default value is 6. It has not be given to the constructor.}

\item{\code{dimension} }{ An integer that declares the dimension of the convex body. It has not be given to the constructor.}

\item{\code{volume} }{ The volume of the convex body, if it is known.}
 }
}
\keyword{internal}

