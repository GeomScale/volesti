#' @title Volesti 2D polytope sampling demo
#' @description Installs required R packages, samples points from a 2D
#' H-polytope using volesti, computes its inner ball and volume, and
#' produces illustrative feasibility and optimization plots.

#set CRAN mirror for non-interactive installs
options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("volesti", "ggplot2", "xfun", "gMOIP", "gridExtra")
for (p in pkgs) {
   if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
   }
}

library(volesti)
library(ggplot2)
library(xfun)
library(gMOIP)
library(gridExtra)

A <- matrix(c(
  1, 0,
 -1, 0,
  0, 1,
  0, -1
), nrow=4, byrow=TRUE)

A_mat <- as.matrix(A)
b <- c(1, 1, 1, 1)
b_vec <- as.vector(b)
objective <- c(0,0)

N <- 10000
P <- volesti::Hpolytope(A = A_mat, b = b_vec)
samples <- volesti::sample_points(P, n = N)
df_sample <- as.data.frame(t(samples))
ball <- inner_ball(P) # centre of the largest inscribed polytope
print(ball)
volume <- volesti::volume(P) # volume of the polytope
print(volume)

ggplot(df_sample, aes(x = V1, y = V2)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  theme_minimal() +
  coord_fixed() +
  labs(title = "2D Projection of Hpolytope Samples", x = "Dim 1", y = "Dim 2")

# https://cran.r-project.org/web/packages/gMOIP/vignettes/polytope_2d.html
p1 <- plotPolytope(A = A,b = b, objective, type = rep("c", ncol(A)), 
    crit = "max", faces = rep("c", ncol(A)), 
    plotFaces = TRUE, plotFeasible = TRUE, 
    plotOptimum = FALSE, labels = NULL
    ) + ggplot2::ggtitle("Feasible region only")

p2 <- plotPolytope(A, b, objective,
   type = rep("c", ncol(A)), crit = "max", faces = rep("c", ncol(A)),
   plotFaces = TRUE, plotFeasible = TRUE, plotOptimum = TRUE, labels = "coord"
) + ggplot2::ggtitle("Solution LP max")

p3 <- plotPolytope(A, b, objective, type = rep("c", ncol(A)),
   crit = "min", faces = rep("c", ncol(A)), plotFaces = TRUE,
   plotFeasible = TRUE, plotOptimum = TRUE, labels = "n"
) + ggplot2::ggtitle("Solution LP min")

p4 <- plotPolytope(A, b, objective,
   type = rep("c", ncol(A)), crit = "max", faces = rep("c", ncol(A)),
   plotFaces = TRUE, plotFeasible = TRUE, plotOptimum = TRUE, labels = "coord"
) + ggplot2::xlab("x") + ggplot2::ylab("y") + ggplot2::ggtitle("Solution (max) with other axis labels")

g <- gridExtra::arrangeGrob(p1, p2, p3, p4, nrow = 2)
ggplot2::ggsave("poly_grid.png", g, width = 12, height = 8)
