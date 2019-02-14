\name{Zonotope}
\alias{Zonotope}
\title{An exposed C++ class to represent zonotopes}

\description{
A zonotope is a convex polytope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
}
\section{Fields}{
\itemize{
\item{\code{G}}{\eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise}
  
\item{\code{type}}{An integer that declares the representation of the polytope. For H-representation the default value is 3. It has not be given to the constructor.}
  
\item{\code{dimension}}{An integer that declares the dimension of the polytope. It has not be given to the constructor.}
}}
