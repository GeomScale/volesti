# volesti benchmarking: Convex Hull Approximation

## What This Demonstrates

volesti excels at sampling **from** known convex bodies. This benchmarking example explores the inverse application: measuring how accurately we can reconstruct a geometric body **from** those samples. 

This example provides a concrete quantitative evaluation of approximation quality using support functions as a proxy for the Hausdorff distance.

1. **Inner approximation via convex hull**: Sample N points from a known body, and theoretically treat their convex hull as an inner approximation.
2. **Support function gap measurement**: For random directions $u$, compute the gap $h_K(u) - \hat{h}_n(u)$ between the true and estimated support functions. The maximum gap tightly approximates the Hausdorff distance.
3. **Convergence rate validation**: Show that the gap decreases as N grows.
4. **Body type auto-detection**: Track how the number of hull vertices grows with N. Polyhedral bodies show logarithmic growth; smooth bodies show polynomial growth.
5. **Holdout miss rate**: Estimate the fraction of the true body's volume roughly not covered by the convex hull.

## Mathematical Background

For a convex body $K \subset \mathbb{R}^d$ and i.i.d. samples $X_1,...,X_n \sim \text{Uniform}(K)$:

- **Support function**: $h_K(u) = \max\{\langle u,x\rangle : x \in K\}$ uniquely characterizes $K$.
- **Hausdorff distance**: $d_H(K, L) = \|h_K - h_L\|_\infty$.
- **Convergence rates**:
  - Polytopes: $d_H(K, \hat{K}_n) = O((\log n / n)^{1/d})$
  - Smooth $C^2$ bodies: $d_H(K, \hat{K}_n) = O((\log n / n)^{2/(d+1)})$
- **Vertex count growth**:
  - Polytopes: $\mathbb{E}[\#\text{vertices}] \sim (\log n)^{(d-1)/2}$
  - Smooth: $\mathbb{E}[\#\text{vertices}] \sim n^{(d-1)/(d+1)}$

## Build and Run

```bash
cd examples/geometric_inference
mkdir build && cd build
cmake ..
make
./geometric_inference_poc
```
