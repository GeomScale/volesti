\docType{class}
\name{Rcpp_Vpolytope}
\alias{Rcpp_Vpolytope-class}
\alias{[,Rcpp_Vpolytope-method}
\alias{[,Rcpp_Vpolytope,ANY,ANY,ANY-method}
\alias{$<-,Rcpp_Vpolytope-method}
\alias{$,Rcpp_Vpolytope-method}
\alias{filepaths<-,Rcpp_Vpolytope-method}
\title{
An \code{Rcpp} class to represent V-polytopes, exposed to \code{R} via modules.
}
\description{
A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which corresponds to its vertices.
}
\details{
\describe{
   \item{\code{V} }{ A \eqn{m\times d} numerical matrix that contains the vertices row-wise}
  
  \item{\code{type} }{ An integer that declares the representation of the polytope. For V-representation the default value is 2.}
  
  \item{\code{dimension} }{ The dimension of the polytope.}
  
  \item{\code{volume} }{ The volume of the polytope, if it is known.}
 }
}
\keyword{internal}
