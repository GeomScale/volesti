\name{Vpolytope}
\alias{Vpolytope}
\title{An \code{R} class to represent V-polytopes.}

\description{
  A V-polytope is defined as the convex hull of \eqn{m} \eqn{d}-dimensional points which corresponds to its vertices.
}
\section{Fields}{
\itemize{
  \item{\code{V} }{ A \eqn{m\times d} numerical matrix that contains the vertices row-wise}
  
  \item{\code{type} }{ An integer that declares the representation of the polytope. For V-representation the default value is 2.}
  
  \item{\code{dimension} }{ The dimension of the polytope.}
  
  \item{\code{volume} }{ The volume of the polytope, if it is known.}
}}
