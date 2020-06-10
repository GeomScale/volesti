\name{Hpolytope}
\alias{Hpolytope}
\title{An \code{R} class to represent H-polytopes.}

\description{
A H-polytope is a convex polytope defined by a set of linear inequalities or equivalently a \eqn{d}-dimensional H-polytope with \eqn{m} facets is defined by a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b, s.t.: \eqn{P=\{x\ |\  Ax\leq b\} }.
}
\section{Fields}{
\itemize{
\item{\code{A} }{ A \eqn{m\times d} numerical matrix A}

\item{\code{b} }{ \eqn{m}-dimensional vector b}

\item{\code{type} }{ An integer that declares the representation of the polytope. For H-representation the default value is 1.}

\item{\code{dimension} }{ The dimension of the polytope.}

\item{\code{volume} }{ The volume of the polytope, if it is known.}
}}
