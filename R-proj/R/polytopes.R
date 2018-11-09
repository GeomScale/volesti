#' A reference class to represent a H-polytope
#' 
#' @description A H-polytope is a convex polytope defined by a set of linear inequalities or equivalently a \eqn{d}-dimensional H-polytope with \eqn{m} facets is defined by a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b, s.t.: \eqn{Ax\leq b}.
#' 
#' @field A \eqn{m\times d} numerical matrix A
#' @field A \eqn{m}-dimensional vector b
HPolytope <- setRefClass("HPolytope", fields = list(A="matrix", b="vector", t="numeric"),
            methods = list(
              get_mat = function() {
                "Return a numerical matrix that describes the polytope in ine format"
                Mat = -A
                vec = b
                d = dim(Mat)[2] + 1
                m = dim(Mat)[1]
                r = rep(0,d)
                r[1] = m
                r[2] = d
                Mat = matrix(cbind(b, Mat), ncol = dim(Mat)[2] + 1)
                Mat = matrix(rbind(r, Mat), ncol = dim(Mat)[2])
                return(Mat)
              }
            ))
#' A reference class to represent a V-polytope or the intersection of two V-polytopes
#' 
#' @description A V-polytope is a convex polytope defined by the set of its vertices.
#' 
#' @field A \eqn{m\times d} numerical matrix V that contains row-wise the \eqn{m} vertices of a V-polytope.
#' @field A \eqn{k\times d} numerical matrix V2 that contains row-wise the \eqn{k} vertices of the second V-polytope that defines the intersection. If this field is NULL then the created object will represent a V-polytope defined by the first field.
VPolytope <- setRefClass("VPolytope", fields = list(V="matrix", V2="matrix", t="numeric"),
            methods = list(
              get_mat = function() {
                "Return a numerical matrix that describes the first polytope in ext format"
                Mat = V
                d = dim(Mat)[2] + 1
                m = dim(Mat)[1]
                b = rep(1, m)
                r = rep(0, d)
                r[1] = m
                r[2] = d
                Mat = matrix(cbind(b, Mat), ncol = dim(Mat)[2] + 1)
                Mat = matrix(rbind(r, Mat), ncol = dim(Mat)[2])
                return(Mat)
              },
              get_mat2 = function() {
                "Return a numerical matrix that describes the second polytope in ine format"
                Mat = V2
                d = dim(Mat)[2] + 1
                m = dim(Mat)[1]
                b = rep(1, m)
                r = rep(0, d)
                r[1] = m
                r[2] = d
                Mat = matrix(cbind(b, Mat), ncol = dim(Mat)[2] + 1)
                Mat = matrix(rbind(r, Mat), ncol = dim(Mat)[2])
                return(Mat)
              }
            ))
#' A class that represent a zonotope
#' 
#' @description A \eqn{d}-dimensional zonotope is a convex polytope defined by the Minkowski sum of \eqn{d}-dimensional segments.
#' 
#' @field A \eqn{m\times d} numerical matrix G that contains row-wise the \eqn{m} segments.
Zonotope <- setRefClass("Zonotope", fields = list(G="matrix", t="numeric"),
            methods = list(
              get_mat = function() {
                "Return a numerical matrix that describes the zonotope in ine format"
                Mat = G
                d = dim(Mat)[2] + 1
                m = dim(Mat)[1]
                b = rep(1, m)
                r = rep(0, d)
                r[1] = m
                r[2] = d
                Mat = matrix(cbind(b, Mat), ncol = dim(Mat)[2] + 1)
                Mat = matrix(rbind(r, Mat), ncol = dim(Mat)[2])
                return(Mat)
              }
            ))