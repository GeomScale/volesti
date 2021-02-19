//Copyright 2018* United States Secretary of Commerce, NIST.
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Credit to https://github.com/usnistgov/ChebTools

#ifndef CHEBTOOLS_H
#define CHEBTOOLS_H

#include "Eigen/Dense"
#include <vector>
#include <queue>

#include <map>


namespace ChebTools{

    // https://proquest.safaribooksonline.com/9780321637413
    // https://web.stanford.edu/class/archive/cs/cs107/cs107.1202/lab1/
    static int midpoint_Knuth(int x, int y) {
        return (x & y) + ((x ^ y) >> 1);
    };


    /**
    * @brief This class stores sets of L matrices (because they are a function only of the degree of the expansion)
    *
    * The L matrix is used to convert from functional values to coefficients, as in \f[ \vec{c} = \mathbf{L}\vec{f} \f]
    */
    class LMatrixLibrary {
    private:
        std::map<std::size_t, Eigen::MatrixXd> matrices;
        void build(std::size_t N) {
            Eigen::MatrixXd L(N + 1, N + 1); ///< Matrix of coefficients
            for (int j = 0; j <= N; ++j) {
                for (int k = j; k <= N; ++k) {
                    double p_j = (j == 0 || j == N) ? 2 : 1;
                    double p_k = (k == 0 || k == N) ? 2 : 1;
                    L(j, k) = 2.0 / (p_j*p_k*N)*cos((j*EIGEN_PI*k) / N);
                    // Exploit symmetry to fill in the symmetric elements in the matrix
                    L(k, j) = L(j, k);
                }
            }
            matrices[N] = L;
        }
    public:
        /// Get the \f$\mathbf{L}\f$ matrix of degree N
        const Eigen::MatrixXd & get(std::size_t N) {
            auto it = matrices.find(N);
            if (it != matrices.end()) {
                return it->second;
            }
            else {
                build(N);
                return matrices.find(N)->second;
            }
        }
    };
    static LMatrixLibrary l_matrix_library;

    /**
    * @brief This class stores sets of U matrices (because they are a function only of the degree of the expansion)
    *
    * The U matrix is used to convert from coefficients to functional values, as in \f[ \vec{f} = \mathbf{U}\vec{c} \f]
    */
    class UMatrixLibrary {
    private:
        std::map<std::size_t, Eigen::MatrixXd> matrices;
        void build(std::size_t N) {
            Eigen::MatrixXd U(N + 1, N + 1); ///< Matrix of coefficients
            for (int j = 0; j <= N; ++j) {
                for (int k = j; k <= N; ++k) {
                    U(j, k) = cos((j*EIGEN_PI*k) / N);
                    // Exploit symmetry to fill in the symmetric elements in the matrix
                    U(k, j) = U(j, k);
                }
            }
            matrices[N] = U;
        }
    public:
        /// Get the \f$\mathbf{U}\f$ matrix of degree N
        const Eigen::MatrixXd & get(std::size_t N) {
            auto it = matrices.find(N);
            if (it != matrices.end()) {
                return it->second;
            }
            else {
                build(N);
                return matrices.find(N)->second;
            }
        }
    };
    static UMatrixLibrary u_matrix_library;

    /**
    For a monotonically increasing vector, find the left index of the interval bracketing the given value
    */
    template<typename VecType>
    int get_increasingleftofval(const VecType& breakpoints, double x, int N) {
        int iL = 0, iR = N - 1, iM;
        while (iR - iL > 1) {
            iM = midpoint_Knuth(iL, iR);
            if (x >= breakpoints[iM]) {
                iL = iM;
            }
            else {
                iR = iM;
            }
        }
        return iL;
    };

    /**
    For a monotonically decreasing vector, find the left index of the interval bracketing the given value
    */
    template<typename VecType>
    int get_decreasingleftofval(const VecType& breakpoints, double x, int N) {
        int iL = 0, iR = N - 1, iM;
        while (iR - iL > 1) {
            iM = midpoint_Knuth(iL, iR);
            if (x <= breakpoints[iM]) {
                iL = iM;
            }
            else {
                iR = iM;
            }
        }
        return iL;
    };

    typedef Eigen::VectorXd vectype;
    
    /// Get the Chebyshev-Lobatto nodes for an expansion of degree \f$N\f$
    const Eigen::VectorXd &get_CLnodes(std::size_t N);

    Eigen::VectorXcd eigenvalues(const Eigen::MatrixXd &A, bool balance);
    Eigen::VectorXd eigenvalues_upperHessenberg(const Eigen::MatrixXd &A, bool balance);

    /**
    * @brief This is the main underlying object that makes all of the code of ChebTools work.
    *
    * This class has accessor methods for getting things from the object, and static factory
    * functions for generating new expansions.  It also has methods for calculating derivatives,
    * roots, etc.
    */
    class ChebyshevExpansion {
    private:
        vectype m_c;
        double m_xmin, m_xmax;

        vectype m_recurrence_buffer;
        vectype m_nodal_value_cache;
        void resize() {
            m_recurrence_buffer.resize(m_c.size());
        }

        //reduce_zeros changes the m_c field so that our companion matrix doesnt have nan values in it
        //all this does is truncate m_c such that there are no trailing zero values
        static Eigen::VectorXd reduce_zeros(const Eigen:: VectorXd &chebCoeffs){
          //these give us a threshold for what coefficients are large enough
          double largeTerm = 1e-15;
          if (chebCoeffs.size()>=1 && std::abs(chebCoeffs(0))>largeTerm){
            largeTerm = chebCoeffs(0);
          }
          //if the second coefficient is larger than the first, then make our tolerance
          //based on the second coefficient, this is useful for functions whose mean value
          //is zero on the interval
          if (chebCoeffs.size()>=2 && std::abs(chebCoeffs(1))>largeTerm){
            largeTerm = chebCoeffs(1);
          }
          double tol = largeTerm*(1e-15);
          int neededSize = static_cast<int>(chebCoeffs.size());
          //loop over m_c backwards, if we run into large enough coefficient, then record the size and break
          for (int i=static_cast<int>(chebCoeffs.size())-1; i>=0; i--){
            if (std::abs(chebCoeffs(i))>tol){
              neededSize = i+1;
              break;
            }
            neededSize--;
          }
          //neededSize gives us the number of coefficients that are nonzero
          //we will resize m_c such that there are essentially no trailing zeros
          return chebCoeffs.head(neededSize);
        }

    public:
        /// Initializer with coefficients, and optionally a range provided
        ChebyshevExpansion(const vectype &c, double xmin = -1, double xmax = 1) : m_c(c), m_xmin(xmin), m_xmax(xmax) { resize(); };
        /// Initializer with coefficients, and optionally a range provided
        ChebyshevExpansion(const std::vector<double> &c, double xmin = -1, double xmax = 1) : m_xmin(xmin), m_xmax(xmax) {
            m_c = Eigen::Map<const Eigen::VectorXd>(&(c[0]), c.size());
            resize();
        };
        /// Move constructor (C++11 only)
        ChebyshevExpansion(const vectype &&c, double xmin = -1, double xmax = 1) : m_c(c), m_xmin(xmin), m_xmax(xmax) { resize(); };

        /// Cache nodal function values
        void cache_nodal_function_values(vectype values) {
            m_nodal_value_cache = values;
        }
        /// Get the minimum value of \f$x\f$ for the expansion
        double xmin() const{ return m_xmin; }
        /// Get the maximum value of \f$x\f$ for the expansion
        double xmax() const{ return m_xmax; }
        /// Go from a value in [xmin,xmax] to a value in [-1,1]
        double scale_x(const double x) const {
            return (2 * x - (m_xmax + m_xmin)) / (m_xmax - m_xmin);
        }
        /// Map from a value in [-1,1] to a value in [xmin,xmax]
        double unscale_x(const double xscaled) const {
            return ((m_xmax - m_xmin)*xscaled + (m_xmax + m_xmin))/2;
        }

        /// Get the vector of coefficients in increasing order
        const vectype &coef() const;

        /// Return the N-th derivative of this expansion, where N must be >= 1
        ChebyshevExpansion deriv(std::size_t Nderiv) const;
        /// Return the indefinite integral of this function
        ChebyshevExpansion integrate(std::size_t Nintegral = 1) const;
        /// Get the Chebyshev-Lobatto nodes in the domain [-1,1]
        Eigen::VectorXd get_nodes_n11();
        /// Get the Chebyshev-Lobatto nodes in the domain [-1,1]; thread-safe const variant
        Eigen::VectorXd get_nodes_n11() const {
            Eigen::Index N = m_c.size() - 1;
            double NN = static_cast<double>(N);
            return (Eigen::VectorXd::LinSpaced(N + 1, 0, NN).array() * EIGEN_PI / N).cos();
        }
        /// Get the Chebyshev-Lobatto nodes in the domain [xmin, xmax]
        Eigen::VectorXd get_nodes_realworld();
        /// Get the Chebyshev-Lobatto nodes in the domain [xmin, xmax]; thread-safe const variant
        Eigen::VectorXd get_nodes_realworld() const {
            return ((m_xmax - m_xmin) * get_nodes_n11().array() + (m_xmax + m_xmin)) * 0.5;
        }

        /// Values of the function at the Chebyshev-Lobatto nodes
        Eigen::VectorXd get_node_function_values() const;
        /// Return true if the function values at the Chebyshev-Lobatto nodes are monotonic with the independent variable
        bool is_monotonic() const;

        // ******************************************************************
        // ***********************      OPERATORS     ***********************
        // ******************************************************************

        /// A ChebyshevExpansion plus another ChebyshevExpansion yields a new ChebyheveExpansion
        ChebyshevExpansion operator+(const ChebyshevExpansion &ce2) const ;
        /** 
        * @brief An inplace addition of two expansions
        * @note The lower degree one is right-padded with zeros to have the same degree as the higher degree one
        * @param donor The other expansion in the summation
        */
        ChebyshevExpansion& operator+=(const ChebyshevExpansion &donor);
        /// Multiplication of an expansion by a constant
        ChebyshevExpansion operator*(double value) const;
        /// Addition of a constant to an expansion
        ChebyshevExpansion operator+(double value) const;
        /// Subtraction of a constant from an expansion
        ChebyshevExpansion operator-(double value) const;
        /// An inplace multiplication of an expansion by a constant
        ChebyshevExpansion& operator*=(double value);
        /// An inplace addition of a constant to an expansion
        ChebyshevExpansion& operator+=(double value);
        /// An inplace subtraction of a constant from an expansion
        ChebyshevExpansion& operator-=(double value);
        /// Unary negation operator
        ChebyshevExpansion operator-() const;
        /// An inplace subtraction of an expansion by another expansion
        ChebyshevExpansion& operator-=(const ChebyshevExpansion &ce2);
        /// An inplace subtraction of an expansion by another expansion
        ChebyshevExpansion operator-(const ChebyshevExpansion& ce2) const;
        /**
         * @brief Multiply two Chebyshev expansions together; thanks to Julia code from Bradley Alpert, NIST
         *
         * Converts padded expansions to nodal functional values, functional values are multiplied together,
         * and then inverse transformation is used to return to coefficients of the product
         * @param ce2 The other expansion
         */
        ChebyshevExpansion operator*(const ChebyshevExpansion &ce2) const;

        /**
         * @brief Divide two expansions by each other.  Right's reciprocal is taken, multiplied by this expansion
         *
         * @param ce2 The other expansion
         */
        ChebyshevExpansion operator/(const ChebyshevExpansion& ce2) const {
            return (*this) * ce2.reciprocal();
        }
        /**
         * @brief Multiply a Chebyshev expansion by its independent variable \f$x\f$
         */
        ChebyshevExpansion times_x() const;

        /** 
         * @brief Multiply a Chebyshev expansion by its independent variable \f$x\f$ in-place
         *
         * This operation is carried out in-place to minimize the amount of memory re-allocation
         * which proved during profiling to be a major source of inefficiency
         */
        ChebyshevExpansion& times_x_inplace();

        ChebyshevExpansion reciprocal() const;

        /// Friend function that allows for pre-multiplication by a constant value
        friend ChebyshevExpansion operator*(double value, const ChebyshevExpansion &ce){
            return ChebyshevExpansion(std::move(ce.coef()*value),ce.m_xmin, ce.m_xmax);
        };
        /// Friend function that allows expansion to be the denominator in division with double
        friend ChebyshevExpansion operator/(double value, const ChebyshevExpansion& ce) {
            return value * ce.reciprocal();
        };
        /// Friend function that allows pre-subtraction of expansion (value-expansion)
        friend ChebyshevExpansion operator-(double value, const ChebyshevExpansion& ce) {
            return -ce+value;
        };
        /// Friend function that allows pre-addition of expansion (value+expansion)
        friend ChebyshevExpansion operator+(double value, const ChebyshevExpansion& ce) {
            return ce + value;
        };
        
        /**
         * @brief Apply a function to the expansion
         *
         * This function first converts the expansion to functional values at the 
         * Chebyshev-Lobatto nodes, applies the function to the nodal values, and then
         * does the inverse transformation to arrive at the coefficients of the expansion
         * after applying the transformation
         */
        ChebyshevExpansion apply(std::function<Eigen::ArrayXd(const Eigen::ArrayXd &)> &f) const;

        // ******************************************************************
        // **********************      EVALUATORS     ***********************
        // ******************************************************************

        /**
        * @brief Do a single input/single output evaluation of the Chebyshev expansion with the inputs scaled in [xmin, xmax]
        * @param x A value scaled in the domain [xmin,xmax]
        */
        double y_recurrence(const double x);
        /**
        * @brief Do a single input/single output evaluation of the Chebyshev expansion with the inputs scaled in [xmin, xmax]
        * @param x A value scaled in the domain [xmin,xmax]
        */
        double y_Clenshaw(const double x) const { return y_Clenshaw_xscaled(scale_x(x));  }
        /**
        * @brief Do a single input/single output evaluation of the Chebyshev expansion with the inputs scaled in [-1,1]
        * @param x A value scaled in the domain [-1,1]
        */
        double y_Clenshaw_xscaled(const double x) const;
        /**
        * @brief Do a vectorized evaluation of the Chebyshev expansion with the inputs scaled in [xmin, xmax]
        * @param x A vectype of values in the domain [xmin,xmax]
        */
        vectype y(const vectype &x) const;
        /**
        * @brief Do a vectorized evaluation of the Chebyshev expansion with the inputs scaled in [xmin, xmax]
        * @param x A value scaled in the domain [xmin,xmax]
        */
        double y(const double x) const{ return y_Clenshaw(x); }
        /**
        * @brief Do a vectorized evaluation of the Chebyshev expansion with the input scaled in the domain [-1,1]
        * @param xscaled A vectype of values scaled to the domain [-1,1] (the domain of the Chebyshev basis functions)
        * @returns y A vectype of values evaluated from the expansion
        *
        * By using vectorizable types like Eigen::MatrixXd, without
        * any additional work, "magical" vectorization is happening
        * under the hood, giving a significant speed improvement. From naive
        * testing, the increase was a factor of about 10x.
        */
        vectype y_recurrence_xscaled(const vectype &xscaled) const ;
        /**
        * @brief Do a vectorized evaluation of the Chebyshev expansion with the input scaled in the domain [-1,1] with Clenshaw's method
        * @param xscaled A vectype of values scaled to the domain [-1,1] (the domain of the Chebyshev basis functions)
        * @returns y A vectype of values evaluated from the expansion
        */
        vectype y_Clenshaw_xscaled(const vectype &xscaled) const ;

        /**
        * @brief Construct and return the companion matrix of the Chebyshev expansion
        * @returns A The companion matrix of the expansion
        *
        * See Boyd, SIAM review, 2013, http://dx.doi.org/10.1137/110838297, Appendix A.2
        */
        Eigen::MatrixXd companion_matrix(const Eigen::VectorXd &coeffs) const ;
        /**
        * @brief Return the real roots of the Chebyshev expansion
        * @param only_in_domain If true, only real roots that are within the domain
        *                       of the expansion will be returned, otherwise all real roots
        *
        * The roots are obtained based on the fact that the eigenvalues of the
        * companion matrix are the roots of the Chebyshev expansion.  Thus
        * this function is relatively slow, because an eigenvalue solve is required,
        * which takes O(n^3) FLOPs.  But it is numerically rather reliable.
        *
        * As the order of the expansion increases, the eigenvalue solver in Eigen becomes
        * progressively less and less able to obtain the roots properly. The eigenvalue
        * solver in numpy tends to be more reliable.
        */
        std::vector<double> real_roots(bool only_in_domain = true) const ;
        /**
        * @brief The second-generation rootfinder of ChebyshevExpansions
        * @param only_in_domain True: only keep roots that are in the domain of the expansion. False: all real roots
        */
        std::vector<double> real_roots2(bool only_in_domain = true) const;

        /**
        * @brief Calculate the value (only one) of x in [xmin, xmax] for which the expansion value is equal to given value
        *
        * Functionally the use is similar to real_roots except that:
        * 1) nodal values are cached
        * 2) only one solution is possible
        *
        Warning: the monotonicity of the expansion is assumed, but not checked
        *
        * @param yval Given value for which value of x is to be obtained
        */
        double monotonic_solvex(double yval);

        /**
        * @brief Subdivide the original interval into a set of subintervals that are linearly spaced
        * @note A vector of ChebyshevExpansions are returned
        * @param Nintervals The number of intervals
        * @param Ndegree The degree of the Chebyshev expansion in each interval
        */
        std::vector<ChebyshevExpansion> subdivide(std::size_t Nintervals, std::size_t Ndegree) const ;

        /**
        * @brief For a vector of ChebyshevExpansions, find all roots in each interval
        * @param segments The vector of ChebyshevExpansions
        * @param only_in_domain True: only keep roots that are in the domain of the expansion. False: all real roots
        */
        static std::vector<double> real_roots_intervals(const std::vector<ChebyshevExpansion> &segments, bool only_in_domain = true);

        /**
        * @brief Time how long (in seconds) it takes to evaluate the roots
        * @param N How many repeats to do (maybe a million?  It's pretty fast for small degrees)
        */
        double real_roots_time(long N);

        /// A DEPRECATED function for approximating the roots (do not use)
        std::vector<double> real_roots_approx(long Npoints);

        // ******************************************************************
        // ***********************      BUILDERS      ***********************
        // ******************************************************************

        /**
        * @brief Given a set of values at the Chebyshev-Lobatto nodes, perhaps obtained from the ChebyshevExpansion::factory function, 
        * get the expansion, using the discrete cosine transform (DCT) approach
        *
        * @param N The degree of the expansion
        * @param f The set of values at the Chebyshev-Lobatto nodes
        * @param xmin The minimum value of x for the expansion
        * @param xmax The maximum value of x for the expansion
        */
        static ChebyshevExpansion factoryf(const std::size_t N, const Eigen::VectorXd &f, const double xmin, const double xmax) ;

        /**
        * @brief Given a set of values at the Chebyshev-Lobatto nodes, build the expansion, using the FFT approach
        *
        * See this clear example: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/23972/versions/22/previews/chebfun/examples/approx/html/ChebfunFFT.html
        *
        * @param N The degree of the expansion
        * @param f The set of values at the Chebyshev-Lobatto nodes
        * @param xmin The minimum value of x for the expansion
        * @param xmax The maximum value of x for the expansion
        */
        static ChebyshevExpansion factoryfFFT(const std::size_t N, const Eigen::VectorXd& f, const double xmin, const double xmax);

        /**
        * @brief Given a callable function, construct the N-th order Chebyshev expansion in [xmin, xmax]
        * @param N The order of the expansion; there will be N+1 coefficients
        * @param func A callable object, taking the x value (in [xmin,xmax]) and returning the y value
        * @param xmin The minimum x value for the fit
        * @param xmax The maximum x value for the fit
        *
        * See Boyd, SIAM review, 2013, http://dx.doi.org/10.1137/110838297, Appendix A.
        */
        template<class double_function>
        static ChebyshevExpansion factory(const std::size_t N, double_function func, const double xmin, const double xmax)
        {
            // Get the precalculated Chebyshev-Lobatto nodes
            const Eigen::VectorXd & x_nodes_n11 = get_CLnodes(N);

            // Step 1&2: Grid points functional values (function evaluated at the
            // extrema of the Chebyshev polynomial of order N - there are N+1 of them)
            Eigen::VectorXd f(N + 1);
            for (int k = 0; k <= N; ++k) {
                // The extrema in [-1,1] scaled to real-world coordinates
                double x_k = ((xmax - xmin)*x_nodes_n11(k) + (xmax + xmin)) / 2.0;
                f(k) = func(x_k);
            }
            return factoryf(N, f, xmin, xmax);
        };

        /// Convert a monomial term in the form \f$x^n\f$ to a Chebyshev expansion
        static ChebyshevExpansion from_powxn(const std::size_t n, const double xmin, const double xmax);

        /** 
        * @brief Convert a polynomial expansion in monomial form to a Chebyshev expansion
        *
        * The monomial expansion is of the form \f$ y = \displaystyle\sum_{i=0}^N c_ix_i\f$
        *
        * This transformation can be carried out analytically.  For convenience we repetitively use
        * calls to ChebyshevExpansion::from_powxn to build up the expansion.  This is probably not
        * the most efficient option, but it works.
        *
        * @param c The vector of coefficients of the monomial expansion in *increasing* degree: 
        * @param xmin The minimum value of \f$x\f$ for the expansion
        * @param xmax The maximum value of \f$x\f$ for the expansion
        */
        template<class vector_type>
        static ChebyshevExpansion from_polynomial(vector_type c, const double xmin, const double xmax) {
            vectype c0(1); c0 << 0;
            ChebyshevExpansion s(c0, xmin, xmax);
            for (std::size_t i = 0; i < static_cast<std::size_t>(c.size()); ++i) {
                s += c(i)*from_powxn(i, xmin, xmax);
            }
            return s;
        }

        static auto dyadic_splitting(const std::size_t N, const std::function<double(double)>& func, const double xmin, const double xmax, 
            const int M, const double tol, const int max_refine_passes = 8, 
            const std::function<void(int, const std::deque<ChebyshevExpansion>&)>&callback = {}) 
        {
            
            // Convenience function to get the M-element norm
            auto get_err = [M](const ChebyshevExpansion& ce) { return ce.coef().tail(M).norm() / ce.coef().head(M).norm(); };
            
            // Start off with the full domain from xmin to xmax
            std::deque<ChebyshevExpansion> expansions;
            expansions.emplace_back(ChebyshevExpansion::factory(N, func, xmin, xmax));

            // Now enter into refinement passes
            for (int refine_pass = 0; refine_pass < max_refine_passes; ++refine_pass) {
                bool all_converged = true;
                // Start at the right and move left because insertions will make the length increase
                for (int iexpansion = static_cast<int>(expansions.size())-1; iexpansion >= 0; --iexpansion) {
                    auto& expan = expansions[iexpansion];
                    auto err = get_err(expan);
                    if (err > tol) {
                        // Splitting is required, do a dyadic split
                        auto xmid = (expan.xmin() + expan.xmax()) / 2;
                        auto newleft = ChebyshevExpansion::factory(N, func, expan.xmin(), xmid);
                        auto newright = ChebyshevExpansion::factory(N, func, xmid, expan.xmax());
                        std::swap(expan, newleft);
                        expansions.insert(expansions.begin() + iexpansion+1, newright);
                        all_converged = false;
                    }
                }
                if (all_converged) { break; }
                if (callback != nullptr) {
                    callback(refine_pass, expansions);
                }
            }
            return expansions;
        }
    };

}; /* namespace ChebTools */
#endif
