\name{Zonotope}
\alias{Zonotope}
\title{An \code{R} class to represent zonotopes.}

\description{
A zonotope is a convex polytope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
}
\section{Fields}{
\itemize{
\item{\code{G} }{ A \eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise}
  
\item{\code{type} }{ An integer that declares the representation of the polytope. For zonotopes the default value is 3.}
  
\item{\code{dimension} }{ The dimension of the polytope.}

\item{\code{volume} }{ The volume of the polytope, if it is known.}
}}
