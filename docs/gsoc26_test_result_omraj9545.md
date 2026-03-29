# GSoC26 Test Results - omraj9545

**GitHub:** [https://github.com/omraj9545](https://github.com/omraj9545)

## Environment
- **OS:** Ubuntu 22.04 (WSL/Docker)
- **Compiler:** `g++ 11.4.0`
- **R:** `4.1.2`

---

## Easy Test: C++ interface + maximum ball + using existing functions

### Build
```bash
cd test/build
cmake ..
make
```

### Run
```bash
./inner_ball_example
```

### Output
```text
radius = 1
center = 0 0
```
**Result:** Computed maximum inner ball for unit square using C++ interface.

---

## Easy Test: R interface + maximum ball

### Run
```r
library(volesti)

P <- gen_cube(3, "H")
ball_vec <- inner_ball(P)


d <- ncol(P@A)

center <- ball_vec[1:d]
radius <- ball_vec[d + 1]

print(center)
print(radius)
```

### Output
```text
[1] 0 0 0
[1] 1
```
**Result:** Maximum ball computed successfully in R interface.

---

## Hard Test: nloptr-based max-ball implementation using R interface

### Implementation
```r
library(nloptr)

max_ball_nloptr <- function(A, b, x0 = NULL,
                            algorithm = "NLOPT_LN_COBYLA",
                            xtol_rel = 1e-8, maxeval = 5000) {
  A <- as.matrix(A)
  b <- as.numeric(b)

  m <- nrow(A)
  n <- ncol(A)
  stopifnot(length(b) == m)

  row_norms <- sqrt(rowSums(A^2))

  if (is.null(x0)) x0 <- c(rep(0, n), 0.1)
  stopifnot(length(x0) == n + 1)

  # scalar objective
  eval_f <- function(z) {
    -z[n + 1]  # minimize -r
  }

  # vector of inequalities g(z) <= 0
  eval_g_ineq <- function(z) {
    x <- z[1:n]
    r <- z[n + 1]
    c(as.vector(A %*% x + r * row_norms - b), -r)
  }

  res <- nloptr(
    x0 = x0,
    eval_f = eval_f,
    eval_g_ineq = eval_g_ineq,
    opts = list(
      algorithm = algorithm,
      xtol_rel = xtol_rel,
      maxeval = maxeval,
      print_level = 0
    )
  )

  z <- res$solution
  list(
    center = z[1:n],
    radius = z[n + 1],
    objective = res$objective,
    status = res$status,
    message = res$message,
    raw = res
  )
}

A <- rbind(
  c( 1, 0),
  c(-1, 0),
  c( 0, 1),
  c( 0,-1)
)
b <- c(1,1,1,1)

out <- max_ball_nloptr(A, b)

out$center
out$radius
out$objective
out$status
out$message
print(out$raw)

cat("center:", signif(out$center, 6), "\n")
cat("radius:", signif(out$radius, 10), "\n")
```

### Results & Output
```text
[1]  6.246424e-11 -6.572057e-10
[1] 1
[1] -1
[1] 4
[1] "NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached."

Call:

nloptr(x0 = x0, eval_f = eval_f, eval_g_ineq = eval_g_ineq, opts = list(algorithm = algorithm,
    xtol_rel = xtol_rel, maxeval = maxeval, print_level = 0))


Minimization using NLopt version 2.10.0

NLopt solver status: 4 ( NLOPT_XTOL_REACHED: Optimization stopped because
xtol_rel or xtol_abs (above) was reached. )

Number of Iterations....: 82
Termination conditions:  xtol_rel: 1e-08        maxeval: 5000
Number of inequality constraints:  5
Number of equality constraints:    0
Optimal value of objective function:  -1.00000000043253
Optimal value of controls: 6.246424e-11 -6.572057e-10 1


center: 6.24642e-11 -6.57206e-10
radius: 1
```