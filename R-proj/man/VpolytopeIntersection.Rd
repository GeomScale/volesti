\name{VpolytopeIntersection}
\alias{VpolytopeIntersection}
\title{An \code{R} class to represent the intersection of two V-polytopes.}

\description{
  An intersection of two V-polytopes, \eqn{P_1}, \eqn{P_2}, is defined by the intersection of the two corresponding convex hulls.
}
\section{Fields}{
\itemize{
\item{\code{V1} }{ The numerical matrix that contains the vertices of \eqn{P_1} row-wise.}
  
\item{\code{V2} }{ The numerical matrix that contains the vertices of \eqn{P_2} row-wise.}
  
\item{\code{type} }{ An integer that declares the representation of the polytope. For this polytope the default value is 4.}
  
\item{\code{dimension} }{ The dimension of the polytope.}

\item{\code{volume} }{ The volume of the polytope, if it is known.}
}}
