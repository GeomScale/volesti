#' @export
HPolytope <- setRefClass("HPolytope", fields = list(A="matrix", b="vector"),
            methods = list(
              get_mat = function() {
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

#' @export
VPolytope <- setRefClass("VPolytope", fields = list(V="matrix"),
            methods = list(
              get_mat = function() {
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
              }
            ))


#' @export
Zonotope <- setRefClass("Zonotope", fields = list(G="matrix"),
            methods = list(
              get_mat = function() {
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