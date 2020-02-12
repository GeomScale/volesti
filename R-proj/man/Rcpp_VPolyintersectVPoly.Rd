\docType{class}
\name{Rcpp_VPolyintersectVPoly}
\alias{Rcpp_VPolyintersectVPoly-class}
\alias{[,Rcpp_VPolyintersectVPoly-method}
\alias{[,Rcpp_VPolyintersectVPoly,ANY,ANY,ANY-method}
\alias{$<-,Rcpp_VPolyintersectVPoly-method}
\alias{$,Rcpp_VPolyintersectVPoly-method}
\alias{filepaths<-,Rcpp_VPolyintersectVPoly-method}
\title{
An \code{Rcpp} class to represent the intersection of two V-polytope, exposed to \code{R} via modules.
}
\description{
An intersection of two V-polytopes, \eqn{P_1}, \eqn{P_2}, is defined by the intersection of the two coresponding convex hulls.
}
\details{
\describe{
\item{\code{V1} }{ The numerical matrix that contains the vertices of \eqn{P_1} row-wise.}
  
\item{\code{V2} }{ The numerical matrix that contains the vertices of \eqn{P_2} row-wise.}
  
\item{\code{type} }{ An integer that declares the representation of the polytope. For these kinf of polytopes the default value is 4.}
  
\item{\code{dimension} }{ An integer that declares the dimension of the polytope.}
 }
}
