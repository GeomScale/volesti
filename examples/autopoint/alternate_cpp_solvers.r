options(repos = c(CRAN = "https://cloud.r-project.org"))

for (p in c("volesti", "Rglpk")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(volesti)
library(Rglpk)

max_inner_ball_rglpk <- function(A, b) {
    n <- ncol(A)          
    m <- nrow(A)          

    row_norms <- sqrt(rowSums(A^2))  
    obj <- c(rep(0, n), -1)
    mat <- cbind(A, row_norms)
    dir <- rep("<=", m)
    bounds <- list(
        lower = list(ind = n + 1L,        val = 0),          # r >= 0
        upper = list(ind = seq_len(n + 1), val = rep(Inf, n + 1))
    )
    result <- Rglpk_solve_LP(
        obj    = obj,
        mat    = mat,
        dir    = dir,
        rhs    = b,
        bounds = bounds,
        types  = rep("C", n + 1),   # all continuous
        max    = FALSE               # minimise -r
    )

    if (result$status != 0) {
        warning("Rglpk couldn't find an optimal solution. P.S.: ", result$status)
    }

    list(
        center = result$solution[seq_len(n)],
        radius = result$solution[n + 1],
        status = result$status
    )
}

run_test <- function(label, A, b, expected_radius = NULL) {
    cat("\n~~~~~~~~~~~~~\n")
    cat("Test:", label, "\n")
    res <- max_inner_ball_rglpk(A, b)
    cat("  Center :", round(res$center, 6), "\n")
    cat("  Radius :", round(res$radius, 6), "\n")
    if (!is.null(expected_radius)) {
        cat("  Expected radius:", expected_radius, "\n")
        cat("  Match  :", isTRUE(all.equal(res$radius, expected_radius, tolerance = 1e-4)), "\n")
    }
    invisible(res)
}

A1 <- matrix(c(1,0, -1,0, 0,1, 0,-1), nrow=4, byrow=TRUE)
b1 <- rep(1, 4)
r1 <- run_test("2D unit hypercube [-1,1]^2", A1, b1, expected_radius = 1)

A2 <- rbind(diag(3), -diag(3))
b2 <- rep(1, 6)
r2 <- run_test("3D unit hypercube [-1,1]^3", A2, b2, expected_radius = 1)

A3 <- matrix(c(-1,0, 0,-1, 1,1), nrow=3, byrow=TRUE)
b3 <- c(0, 0, 1)
r3 <- run_test("2D standard simplex", A3, b3, expected_radius = round(1 / (2 + sqrt(2)), 4))


cat("\n~~~~~~~~~~~~~\n")
cat("Cross-validation with volesti::inner_ball (2D cube)\n")
P  <- volesti::Hpolytope(A = A1, b = b1)
vb <- inner_ball(P)
cat("  volesti center:", round(vb[1:2], 6), "  radius:", round(vb[3], 6), "\n")
cat("  Rglpk  center:", round(r1$center, 6), "  radius:", round(r1$radius, 6), "\n")
cat("  Radii match   :", isTRUE(all.equal(vb[3], r1$radius, tolerance = 1e-4)), "\n")


cat("\n~~~~~~~~~~~~~\n")
cat("Test: 10D hypercube (stress test)\n")
A5 <- rbind(diag(10), -diag(10))
b5 <- rep(1, 20)
r5 <- run_test("10D unit hypercube", A5, b5, expected_radius = 1)