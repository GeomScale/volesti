\name{Vpolytope}
\alias{Vpolytope}
\title{An exposed C++ class to represent V-polytopes}

\description{
  A V-polytope is a convex polytope defined by the set of its vertices.
}
\section{Fields}{
\itemize{
  \item{\code{V} }{A \eqn{m\times d} numerical matrix that contains the vertices row-wise}
  
  \item{\code{type} }{An integer that declares the representation of the polytope. For H-representation the default value is 2. It has not be given to the constructor.}
  
  \item{\code{dimension} }{An integer that declares the dimension of the polytope. It has not be given to the constructor.}
}}
