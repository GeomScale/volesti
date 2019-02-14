\name{IntVP}
\alias{IntVP}
\title{An exposed C++ class to represent the intersection of two V-polytopes}

\description{
  An intersection of two V-polytopes is defined by the intersection of the two coresponding convex hulls.
}
\section{Fields}{
\itemize{
\item{\code{V1}}{The numerical matrix that contains the vertices of \eqn{P_1} row-wise.}
  
\item{\code{V2}}{The numerical matrix that contains the vertices of \eqn{P_2} row-wise.}
  
\item{\code{type}}{An integer that declares the representation of the polytope. For H-representation the default value is 1. It has not be given to the constructor.}
  
\item{\code{dimension}}{An integer that declares the dimension of the polytope. It has not be given to the constructor.}
}}
