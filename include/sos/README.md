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

* Dynamically switch between float -> double -> long double -> multiprecision when Matrices become ill-conditioned.
* Parameter tuning
* Choose Method for QR Decompositions dynamically. Benchmarks can be found [here](https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html).
* Higher order corrector steps
* Make sure that no heap memory is allocated dynamically.
* Add Weighted Sum-of-Squares (WSOS) support
* Add multivariate support
* Implement Lagrange polynomials more efficiently (or find C++ library with sufficient precision.)
* Speed up with Intel MKL (considering it is not open source we might not use it !?) or other high-performance library
* Benchmarking with alfonso, SOSTOOLS, MOSEK, SeDuMi 
* Implementation/inspiration from  [4]. In particular interesting are:
    * Positive definite rescaling matrix for added stability and monotonicity 
    * Combining Predictor- and Corrector steps.
* Higher order corrector steps as in [1]
* Write tests
* Include ARPACK/LAPACK. Necessary for sparse Matrix computation. In particular the Polynomial Envelope problem is 
sparse (2 nonzeros per row) which one should be able to exploit.
* Test different central path neighborhoods, e.g. measuring \psi in \infty norm. 
* Interfaces for R (and Python)
    
#### Current issues

* Using double for polynomials of degree > 100 is unstable. On the other hand, high precision 
numbers immensely slow down the IPM
* Generation of Chebyshev Points / Interpolant Basis dominates the runtime.
* High degree polynomial (beginning at around degree 25) lead to oscillation near the boundary. But this might be an artifact from these polynomials as opposed to a bug or lack in precision.

#### Weekly Progress

* Week 1: Implementation of [1,2]. 
* Week 2: Testing implementation with LP and SDP barrier, Debugging, Refactoring.
* Week 3: Implementation of Interpolant-Barrier, Focus on speeding up implementation with dedicated libraries.
* Week 4: Implementation of Toy Problem: Polynomial Envelope.
* Week 5: Discovered [4], Understanding paper and draw ideas to implement more efficiently.
* Week 6: Debugging, Refactoring Code, Speeding up initialisation of "Polynomial Envelope" Problem
 (non-trivial: properties of Chebyshev polynomials, Clenshaw-Curtis algorithm)
 * Week 7: Runtime speedups, more efficient (including experimentally) implementation.
 * Week 8: Further speedup, Big code refactoring, Add Polynomial minimization, Create Benchmarks with alfonso.

#### References

[1] A. Skajaa and Y. Ye, [A homogeneous interior-point algorithm for nonsymmetric convex conic optimization](https://web.stanford.edu/~yyye/nonsymmhsdimp.pdf), Mathematical Programming Ser. A, 150 (2015), pp. 391-422. 

[2] D. Papp and S. Yildiz, On “A homogeneous interior-point algorithm for nonsymmetric convex conic optimization”. https://arxiv.org/abs/1712.00492

[3] D. Papp and S. Yildiz, [Sum-of-squares optimization without semidefinite programming](https://arxiv.org/abs/1712.01792). SIAM Journal on Optimization 29(1), 2019, pp. 822-851. 

[4] R. Badenbroek and J. Dahl, [An Algorithm for Nonsymmetric Conic Optimization Inspired by MOSEK](https://arxiv.org/pdf/2003.01546.pdf). https://arxiv.org/pdf/2003.01546.pdf 

## Advice on Implementation (Elias Talk Friday 17th July)

* static polymorphism
* set up experiment to compare runtimes with MATLAB(alfonso) (univariate or small bivariate):
    * Script for profiling
* FFT to speed up things?
* Gaussian quadrature vs Clenshaw Curtis weight.

* Expansion:
    * DSOS 
* Newton Polytop ("circuits polynomial")
* SONC polynomials: there is no Barrier function. ()
* Hyperbolic polynomials (Renegars algorithm)
* Power cone (Primal and Dual have Barrier functions)

Conclusion Zoom Call 24/07:
* Share report pdf file
* Check if LAPACK would speed up runtimes
* polish PR Request:
	* Reduce external files, prompt the user to provide it
