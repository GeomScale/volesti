### Sum of Squares optimization

This subproject implements the algorithm(s) in [1,2,3]. Example of usage can be found [here](../../doc/cpp_interface.md). 
#### Currently available:

* The generic implementation of the algorithm in [1] with barrier methods for the following cones:
    * Non-negative orthant
    * SDP cone
    * Dual of SOS cone with following bases:
        * monomial
        * interpolant
* Generation of Lagrange polynomials for Chebyshev points of the second kind
* Tool that approximates the polynomial lower envelope of any (stable to degree approx. 100) set of univariate polynomials 
on the interval [-1,1].
* Visualisation via Matplotlib plot.
* Support for several floating-point precision data types for the IPM 
(Primitive Floating Point and boost:multiprecision floats)

#### Next steps:

* Parameter tuning
* Higher order corrector steps
* Add Weighted Sum-of-Squares (WSOS) support
* Add multivariate support
* Implement Lagrange polynomials more efficiently (or find C++ library with sufficient precision.)
* Speed up with Intel MKL (considering it is not open source we miht not use it !?) or other high-performance library
* Implementation/inspiration from  [4]. In particular interesting are:
* Writing tests
* Interfaces for R (and Python)
    
#### Current issues

* Using double for polynomials of degree > 100 is unstable. On the other hand, high precision 
numbers immensely slow down the IPM
* Generation of Chebyshev Points / Interpolant Basis dominates the runtime.
* High degree polynomial (beginning at around degree 25) lead to oscillation near the boundary. But this might be an artifact from these polynomials as opposed to a bug or lack in precision.

#### References

[1] A. Skajaa and Y. Ye, [A homogeneous interior-point algorithm for nonsymmetric convex conic optimization](https://web.stanford.edu/~yyye/nonsymmhsdimp.pdf), Mathematical Programming Ser. A, 150 (2015), pp. 391-422. 

[2] D. Papp and S. Yildiz, On “A homogeneous interior-point algorithm for nonsymmetric convex conic optimization”. https://arxiv.org/abs/1712.00492

[3] D. Papp and S. Yildiz, [Sum-of-squares optimization without semidefinite programming](https://arxiv.org/abs/1712.01792). SIAM Journal on Optimization 29(1), 2019, pp. 822-851. 

[4] R. Badenbroek and J. Dahl, [An Algorithm for Nonsymmetric Conic Optimization Inspired by MOSEK](https://arxiv.org/pdf/2003.01546.pdf). https://arxiv.org/pdf/2003.01546.pdf 

