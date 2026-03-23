#' @title Max-inscribed-ball NLP demo with nloptr
#' @description Formulates and solves the max-inscribed-ball problem for
#' a 2D H-polytope using nloptr by maximizing the radius under linear
#' inequality constraints augmented by row norms.

#set CRAN mirror for non-interactive installs
options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("volesti", "ggplot2", "xfun", "gMOIP", "gridExtra")
for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p)
    }
}

library(ggplot2)
library(nloptr)

A <- matrix(c(
    1, 0,
    -1, 0,
    0, 1,
    0, -1
), nrow=4, byrow=TRUE)

A_mat <- as.matrix(A)
b <- c(1, 1, 1, 1)
eval_f <- function(z) {
    r <- z[length(z)]
    return(-r)  # maximize r = minimize -r
}

eval_grad_f <- function(z) {
    grad <- rep(0, length(z))
    grad[length(z)] <- -1
    return(grad)
}
## constraint: ||ai||^T*x + ai*r - bi <= 0
eval_g_ineq <- function(z, A, b) {
    x <- z[-length(z)]
    r <- z[length(z)]

    constraints <- numeric(nrow(A))
    norms <- sqrt(rowSums(A^2))
    for (i in 1:nrow(A)) {
        ai <- A[i, ]
        constraints[i] <- sum(ai * x) + norms[i] * r - b[i]
    }

    return(constraints)
}

eval_jac_g_ineq <- function(z, A, b) {
    n <- ncol(A)
    m <- nrow(A)
    jac <- matrix(0, nrow = m, ncol = n + 1)
    for (i in 1:m) {
        ai <- A[i, ]
        jac[i, 1:n] <- ai
        jac[i, n + 1] <- sqrt(sum(ai^2))
    }
    return(jac)
}
eps <- 1e-6
# initial vector: n variables for x plus 1 for r
z0 <- c(rep(0, ncol(A)), eps)
cat('length(z0)=', length(z0), '\n')
result <- nloptr(
    x0 = z0,
    eval_f = eval_f,
    eval_grad_f = eval_grad_f,
    eval_g_ineq = function(z) eval_g_ineq(z, A, b),
    eval_jac_g_ineq = function(z) eval_jac_g_ineq(z, A, b),
    opts = list(
        algorithm = "NLOPT_LD_MMA",
        xtol_rel = 1e-6
    )
)

solution <- result$solution
center <- solution[1:ncol(A)]
radius <- solution[ncol(A) + 1]
