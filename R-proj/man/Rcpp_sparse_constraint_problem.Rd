\docType{class}
\name{Rcpp_sparse_constraint_problem}
\alias{Rcpp_sparse_constraint_problem-class}
\alias{[,Rcpp_sparse_constraint_problem-method}
\alias{[,Rcpp_sparse_constraint_problem,ANY,ANY,ANY-method}
\alias{$<-,Rcpp_sparse_constraint_problem-method}
\alias{$,Rcpp_sparse_constraint_problem-method}
\alias{filepaths<-,Rcpp_sparse_constraint_problem-method}
\title{
An \code{Rcpp} class to represent sparse_constraint_problems, exposed to \code{R} via modules.
}
\description{
A constraint problem is defined by a set of linear inequalities and equalities or equivalently a \eqn{d}-dimensional constraint problem is defined by a \eqn{mineq\times d} matrix Aineq and a \eqn{mineq}-dimensional vector bineq, s.t.: \eqn{Aineqx\leq bineq}, a \eqn{meq\times d} matrix Aeq and a \eqn{meq}-dimensional vector beq, s.t.: \eqn{Aeqx= beq} and two \eqn{d} vectors lb, ub such that \eqn{lb\leq x \leq ub}.}
\details{
\describe{
\item{\code{Aineq} }{\eqn{mineq\times d} sparse matrix Aineq}

\item{\code{bineq} }{\eqn{mineq}-dimensional vector bineq}

\item{\code{Aeq} }{\eqn{meq\times d} sparse matrix Aeq}

\item{\code{beq} }{\eqn{meq}-dimensional vector beq}

\item{\code{lb} }{\eqn{d}-dimensional vector bineq}

\item{\code{ub} }{\eqn{d}-dimensional vector bineq}

\item{\code{type} }{An integer that declares the representation of the polytope. For sparse_constraint_problem the default value is 5.}

\item{\code{dimension} }{The dimension of the polytope.}

 }
}
\keyword{internal}
