\docType{class}
\name{Rcpp_Zonotope}
\alias{Rcpp_Zonotope-class}
\alias{[,Rcpp_Zonotope-method}
\alias{[,Rcpp_Zonotope,ANY,ANY,ANY-method}
\alias{$<-,Rcpp_Zonotope-method}
\alias{$,Rcpp_Zonotope-method}
\alias{filepaths<-,Rcpp_Zonotope-method}
\title{
An \code{Rcpp} class to represent zonotopes, exposed to \code{R} via modules.
}
\description{
A zonotope is a convex polytope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
}
\details{
\describe{
   \item{\code{G} }{ A \eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise}
  
\item{\code{type} }{ An integer that declares the representation of the polytope. For zonotopes the default value is 3.}
  
\item{\code{dimension} }{ The dimension of the polytope.}

\item{\code{volume} }{ The volume of the polytope, if it is known.}
 }
}
\keyword{internal}
