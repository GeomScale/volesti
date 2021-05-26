//Copyright 2018* United States Secretary of Commerce, NIST.
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Credit to https://github.com/usnistgov/ChebTools

#include "ChebTools/ChebTools.h"
#include "Eigen/Dense"
#include <unsupported/Eigen/FFT>

#include <algorithm>
#include <functional>

#include <chrono>
#include <iostream>
#include <map>
#include <limits>

#ifndef DBL_EPSILON
#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

#ifndef DBL_MAX
#define DBL_MAX std::numeric_limits<double>::max()
#endif

template <class T> T POW2(T x){ return x*x; }
template <class T> bool inbetween(T bound1, T bound2, T val){ return val > std::min(bound1,bound2) && val < std::max(bound1,bound2); }

/// See python code in https://en.wikipedia.org/wiki/Binomial_coefficient#Binomial_coefficient_in_programming_languages
/// This is a direct translation of that code to C++
double binomialCoefficient(const double n, const double k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    double _k = std::min(k, n - k); //# take advantage of symmetry
    double c = 1;
    for (double i = 0; i < _k; ++i) {
        c *= (n - i) / (i + 1);
    }
    return c;
}

namespace ChebTools {

    inline bool ValidNumber(double x){
        // Idea from http://www.johndcook.com/IEEE_exceptions_in_cpp.html
        return (x <= DBL_MAX && x >= -DBL_MAX);
    };

    void balance_matrix(const Eigen::MatrixXd &A, Eigen::MatrixXd &Aprime, Eigen::MatrixXd &D) {
        // https://arxiv.org/pdf/1401.5766.pdf (Algorithm #3)
        const int p = 2;
        double beta = 2; // Radix base (2?)
        int iter = 0;
        Aprime = A;
        D = Eigen::MatrixXd::Identity(A.rows(), A.cols());
        bool converged = false;
        do {
            converged = true;
            for (Eigen::Index i = 0; i < A.rows(); ++i) {
                double c = Aprime.col(i).lpNorm<p>();
                double r = Aprime.row(i).lpNorm<p>();
                double s = pow(c, p) + pow(r, p);
                double f = 1;
                //if (!ValidNumber(c)){
                //    std::cout << A << std::endl;
                //    throw std::range_error("c is not a valid number in balance_matrix"); }
                //if (!ValidNumber(r)) { throw std::range_error("r is not a valid number in balance_matrix"); }
                while (c < r/beta) {
                    c *= beta;
                    r /= beta;
                    f *= beta;
                }
                while (c >= r*beta) {
                    c /= beta;
                    r *= beta;
                    f /= beta;
                }
                if (pow(c, p) + pow(r, p) < 0.95*s) {
                    converged = false;
                    D(i, i) *= f;
                    Aprime.col(i) *= f;
                    Aprime.row(i) /= f;
                }
            }
            iter++;
            if (iter > 50) {
                break;
            }
        } while (!converged);
    }

    /**
    * @brief A library that stores the Chebyshev-Lobatto nodes in the domain [-1,1]
    * @note The Chebyshev-Lobatto nodes are a function of degree, but not of the coefficients
    */
    class ChebyshevLobattoNodesLibrary {
    private:
        std::map<std::size_t, Eigen::VectorXd> vectors;
        void build(std::size_t N) {
            double NN = static_cast<double>(N); // just a cast
            vectors[N] = (Eigen::VectorXd::LinSpaced(N + 1, 0, NN).array()*EIGEN_PI / N).cos();
        }
    public:
        /// Get the Chebyshev-Lobatto nodes for expansion of degree \f$N\f$
        const Eigen::VectorXd & get(std::size_t N) {
            auto it = vectors.find(N);
            if (it != vectors.end()) {
                return it->second;
            }
            else {
                build(N);
                return vectors.find(N)->second;
            }
        }
    };
    static ChebyshevLobattoNodesLibrary CLnodes_library;
    const Eigen::VectorXd &get_CLnodes(std::size_t N){
        return CLnodes_library.get(N);
    }

    // DEPRECATED: but held around for future reference
    //class ChebyshevRootsLibrary {
    //private:
    //    std::map<std::size_t, Eigen::VectorXd> vectors;
    //    void build(std::size_t N) {
    //        double NN = static_cast<double>(N); // just a cast
    //        vectors[N] = ((Eigen::VectorXd::LinSpaced(N, 0, NN - 1).array() + 0.5)*EIGEN_PI / NN).cos();
    //    }
    //public:
    //    const Eigen::VectorXd & get(std::size_t N) {
    //        auto it = vectors.find(N);
    //        if (it != vectors.end()) {
    //            return it->second;
    //        }
    //        else {
    //            build(N);
    //            return vectors.find(N)->second;
    //        }
    //    }
    //};
    //static ChebyshevRootsLibrary roots_library;

    // From CoolProp
    template<class T> bool is_in_closed_range(T x1, T x2, T x) { return (x >= std::min(x1, x2) && x <= std::max(x1, x2)); };

    ChebyshevExpansion ChebyshevExpansion::operator+(const ChebyshevExpansion &ce2) const {
        if (m_c.size() == ce2.coef().size()) {
            // Both are the same size, nothing creative to do, just add the coefficients
            return ChebyshevExpansion(std::move(ce2.coef() + m_c), m_xmin, m_xmax);
        }
        else{
            if (m_c.size() > ce2.coef().size()) {
                Eigen::VectorXd c(m_c.size()); c.setZero(); c.head(ce2.coef().size()) = ce2.coef();
                return ChebyshevExpansion(c+m_c, m_xmin, m_xmax);
            }
            else {
                std::size_t n = ce2.coef().size();
                Eigen::VectorXd c(n); c.setZero(); c.head(m_c.size()) = m_c;
                return ChebyshevExpansion(c + ce2.coef(), m_xmin, m_xmax);
            }
        }
    };
    ChebyshevExpansion& ChebyshevExpansion::operator+=(const ChebyshevExpansion &donor) {
        std::size_t Ndonor = donor.coef().size(), N1 = m_c.size();
        std::size_t Nmin = std::min(N1, Ndonor), Nmax = std::max(N1, Ndonor);
        // The first Nmin terms overlap between the two vectors
        m_c.head(Nmin) += donor.coef().head(Nmin);
        // If the donor vector is longer than the current vector, resizing is needed
        if (Ndonor > N1) {
            // Resize but leave values as they were
            m_c.conservativeResize(Ndonor);
            // Copy the last Nmax-Nmin values from the donor
            m_c.tail(Nmax - Nmin) = donor.coef().tail(Nmax - Nmin);
        }
        return *this;
    }
    ChebyshevExpansion& ChebyshevExpansion::operator-=(const ChebyshevExpansion& donor) {
        std::size_t Ndonor = donor.coef().size(), N1 = m_c.size();
        std::size_t Nmin = std::min(N1, Ndonor), Nmax = std::max(N1, Ndonor);
        // The first Nmin terms overlap between the two vectors
        m_c.head(Nmin) -= donor.coef().head(Nmin);
        // If the donor vector is longer than the current vector, resizing is needed
        if (Ndonor > N1) {
            // Resize but leave values as they were
            m_c.conservativeResize(Ndonor);
            // Copy the last Nmax-Nmin values from the donor
            m_c.tail(Nmax - Nmin) = -donor.coef().tail(Nmax - Nmin);
        }
        return *this;
    }
    ChebyshevExpansion ChebyshevExpansion::operator-(const ChebyshevExpansion& ce2) const {
        if (m_c.size() == ce2.coef().size()) {
            // Both are the same size, nothing creative to do, just subtract the coefficients
            return ChebyshevExpansion(std::move(ce2.coef() - m_c), m_xmin, m_xmax);
        }
        else {
            if (m_c.size() > ce2.coef().size()) {
                Eigen::VectorXd c(m_c.size()); c.setZero(); c.head(ce2.coef().size()) = ce2.coef();
                return ChebyshevExpansion(m_c - c, m_xmin, m_xmax);
            }
            else {
                std::size_t n = ce2.coef().size();
                Eigen::VectorXd c(n); c.setZero(); c.head(m_c.size()) = m_c;
                return ChebyshevExpansion(c - ce2.coef(), m_xmin, m_xmax);
            }
        }
    };
    ChebyshevExpansion ChebyshevExpansion::operator*(double value) const {
        return ChebyshevExpansion(m_c*value, m_xmin, m_xmax);
    }
    ChebyshevExpansion ChebyshevExpansion::operator+(double value) const {
        Eigen::VectorXd c = m_c;
        c(0) += value;
        return ChebyshevExpansion(c, m_xmin, m_xmax);
    }
    ChebyshevExpansion ChebyshevExpansion::operator-(double value) const {
        Eigen::VectorXd c = m_c;
        c(0) -= value;
        return ChebyshevExpansion(c, m_xmin, m_xmax);
    }
    ChebyshevExpansion ChebyshevExpansion::operator-() const{
        Eigen::VectorXd c = m_c;
        return ChebyshevExpansion(-c, m_xmin, m_xmax);
    }
    ChebyshevExpansion& ChebyshevExpansion::operator*=(double value) {
        m_c *= value;
        return *this;
    }
    ChebyshevExpansion& ChebyshevExpansion::operator+=(double value) {
        m_c(0) += value;
        return *this;
    }
    ChebyshevExpansion& ChebyshevExpansion::operator-=(double value) {
        m_c(0) -= value;
        return *this;
    }
    ChebyshevExpansion ChebyshevExpansion::operator*(const ChebyshevExpansion &ce2) const {

        std::size_t order1 = this->m_c.size()-1,
                    order2 = ce2.coef().size()-1;
        // The order of the product is the sum of the orders of the two expansions
        std::size_t Norder_product = order1 + order2;

        // Create padded vectors, and copy into them the coefficients from this instance
        // and that of the donor
        Eigen::VectorXd a = Eigen::VectorXd::Zero(Norder_product+1),
                        b = Eigen::VectorXd::Zero(Norder_product+1);
        a.head(order1+1) = this->m_c; b.head(order2+1) = ce2.coef();

        // Get the matrices U and V from the libraries
        const Eigen::MatrixXd &U = u_matrix_library.get(Norder_product);
        const Eigen::MatrixXd &V = l_matrix_library.get(Norder_product);

        // Carry out the calculation of the final coefficients
        // U*a is the functional values at the Chebyshev-Lobatto nodes for the first expansion
        // U*b is the functional values at the Chebyshev-Lobatto nodes for the second expansion
        // The functional values are multiplied together in an element-wise sense - this is why both products are turned into arrays
        // The pre-multiplication by V takes us back to coefficients
        return ChebyshevExpansion(V*((U*a).array() * (U*b).array()).matrix(), m_xmin, m_xmax);
    };
    ChebyshevExpansion ChebyshevExpansion::times_x() const {
        // First we treat the of chi*A multiplication in the domain [-1,1]
        Eigen::Index N = m_c.size()-1; // N is the order of A
        Eigen::VectorXd cc(N+2); // Order of x*A is one higher than that of A
        if (N > 1) {
            cc(0) = m_c(1)/2.0;
        }
        if (N > 2) {
            cc(1) = m_c(0) + m_c(2)/2.0;
        }
        for (Eigen::Index i = 2; i < cc.size(); ++i) {
            cc(i) = (i+1 <= N) ? 0.5*(m_c(i-1) + m_c(i+1)) : 0.5*(m_c(i - 1));
        }
        // Scale the values into the real world, which is given by
        // C_scaled = (b-a)/2*(chi*A) + ((b+a)/2)*A
        // where the coefficients in the second term need to be padded with a zero to have
        // the same order as the product of x*A
        Eigen::VectorXd c_padded(N+2); c_padded << m_c, 0;
        Eigen::VectorXd coefs = (((m_xmax - m_xmin)/2.0)*cc).array() + (m_xmax + m_xmin)/2.0*c_padded.array();
        return ChebyshevExpansion(coefs, m_xmin, m_xmax);
    };
    ChebyshevExpansion& ChebyshevExpansion::times_x_inplace() {
        Eigen::Index N = m_c.size() - 1; // N is the order of A
        double diff = ((m_xmax - m_xmin) / 2.0), plus = (m_xmax + m_xmin) / 2.0;
        double cim1old = 0, ciold = 0;
        m_c.conservativeResize(N+2);
        m_c(N+1) = 0.0; // Fill the last entry with a zero
        if (N > 1) {
            // 0-th element
            cim1old = m_c(0); // store the current 0-th element in temporary variable
            m_c(0) = diff*(0.5*m_c(1)) + plus*m_c(0);
        }
        if (N > 2) {
            // 1-th element
            ciold = m_c(1); // store the current 1-th element in temporary variable
            m_c(1) = diff*(cim1old + 0.5*m_c(2)) + plus*m_c(1);
            cim1old = ciold;
        }
        for (Eigen::Index i = 2; i <= N-1; ++i) {
            ciold = m_c(i); // store the current i-th element in temporary variable
            m_c(i) = diff*(0.5*(cim1old + m_c(i + 1)))+plus*m_c(i);
            cim1old = ciold;
        }
        for (Eigen::Index i = N; i <= N + 1; ++i) {
            ciold = m_c(i); // store the current i-th element in temporary variable
            m_c(i) = diff*(0.5*cim1old) + plus*m_c(i);
            cim1old = ciold;
        }
        return *this;
    };
    ChebyshevExpansion ChebyshevExpansion::reciprocal() const{
        // 1. Transform Chebyshev-Lobatto node function values by the function f(y) -> 1/y
        // 2. Go backwards to coefficients from node values c2 = V/y
        const auto Ndegree = m_c.size() - 1;
        const Eigen::MatrixXd& V = l_matrix_library.get(Ndegree);
        // Values at the nodes in the x range of [-1, 1]
        Eigen::VectorXd c = V*(1.0/get_node_function_values().array()).matrix();
        
        return ChebyshevExpansion(c, xmin(), xmax());
    }
    ChebyshevExpansion ChebyshevExpansion::apply(std::function<Eigen::ArrayXd(const Eigen::ArrayXd &)> &f) const{
        // 1. Transform Chebyshev-Lobatto node function values by the function f(y) -> y2
        // 2. Go backwards to coefficients from node values c2 = V*y2
        const auto Ndegree = m_c.size()-1;
        const Eigen::MatrixXd &V = l_matrix_library.get(Ndegree);
        return ChebyshevExpansion(V*f(get_node_function_values()).matrix(), xmin(), xmax());
    }
    bool ChebyshevExpansion::is_monotonic() const {
        auto yvals = get_node_function_values();
        auto N = yvals.size();
        Eigen::ArrayXd diff = yvals.tail(N - 1) - yvals.head(N - 1);
        return (diff < 0.0).all() || (diff > 0.0).all();
    }

    const vectype &ChebyshevExpansion::coef() const {
        return m_c;
    };
    /**
    * @brief Do a single input/single output evaluation of the Chebyshev expansion with the inputs scaled in [xmin, xmax]
    * @param x A value scaled in the domain [xmin,xmax]
    */
    double ChebyshevExpansion::y_recurrence(const double x) {
        // Use the recurrence relationships to evaluate the Chebyshev expansion
        std::size_t Norder = m_c.size() - 1;
        // Scale x linearly into the domain [-1, 1]
        double xscaled = (2 * x - (m_xmax + m_xmin)) / (m_xmax - m_xmin);
        // Short circuit if not using recursive solution
        if (Norder == 0){ return m_c[0]; }
        if (Norder == 1) { return m_c[0] + m_c[1]*xscaled; }

        vectype &o = m_recurrence_buffer;
        o(0) = 1;
        o(1) = xscaled;
        for (int n = 1; n < Norder; ++n) {
            o(n + 1) = 2 * xscaled*o(n) - o(n - 1);
        }
        return m_c.dot(o);
    }
    double ChebyshevExpansion::y_Clenshaw_xscaled(const double xscaled) const {
        // See https://en.wikipedia.org/wiki/Clenshaw_algorithm#Special_case_for_Chebyshev_series
        std::size_t Norder = m_c.size() - 1;
        double u_k = 0, u_kp1 = m_c[Norder], u_kp2 = 0;
        int k = 0;
        for (k = static_cast<int>(Norder) - 1; k >= 1; --k) {
            // Do the recurrent calculation
            u_k = 2.0 * xscaled * u_kp1 - u_kp2 + m_c(k);
            // Update the values
            u_kp2 = u_kp1; u_kp1 = u_k;
        }
        return m_c(0) + xscaled * u_kp1 - u_kp2;
    }
    /**
    * @brief Do a vectorized evaluation of the Chebyshev expansion with the inputs scaled in [xmin, xmax]
    * @param x A vectype of values in the domain [xmin,xmax]
    */
    vectype ChebyshevExpansion::y(const vectype &x) const {
        // Scale x linearly into the domain [-1, 1]
        const vectype xscaled = (2 * x.array() - (m_xmax + m_xmin)) / (m_xmax - m_xmin);
        // Then call the function that takes the scaled x values
        return y_recurrence_xscaled(xscaled);
    }
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
    vectype ChebyshevExpansion::y_recurrence_xscaled(const vectype &xscaled) const {
        const std::size_t Norder = m_c.size() - 1;

        Eigen::MatrixXd A(xscaled.size(), Norder + 1);

        if (Norder == 0) { return m_c[0]*Eigen::MatrixXd::Ones(A.rows(), A.cols()); }
        if (Norder == 1) { return m_c[0] + m_c[1]*xscaled.array(); }

        // Use the recurrence relationships to evaluate the Chebyshev expansion
        // In this case we do column-wise evaluations of the recurrence rule
        A.col(0).fill(1);
        A.col(1) = xscaled;
        for (int n = 1; n < Norder; ++n) {
            A.col(n + 1).array() = 2 * xscaled.array()*A.col(n).array() - A.col(n - 1).array();
        }
        // In this form, the matrix-vector product will yield the y values
        return A*m_c;
    }
    vectype ChebyshevExpansion::y_Clenshaw_xscaled(const vectype &xscaled) const {
        const std::size_t Norder = m_c.size() - 1;
        vectype u_k, u_kp1(xscaled.size()), u_kp2(xscaled.size());
        u_kp1.fill(m_c[Norder]); u_kp2.fill(0);
        for (int k = static_cast<int>(Norder) - 1; k >= 1; --k) {
            u_k = 2 * xscaled.array()*u_kp1.array() - u_kp2.array() + m_c(k);
            // Update summation values for all but the last step
            if (k > 1) {
                u_kp2 = u_kp1; u_kp1 = u_k;
            }
        }
        return xscaled.array()*u_k.array() - u_kp1.array() + m_c(0);
    }

    Eigen::MatrixXd ChebyshevExpansion::companion_matrix(const Eigen::VectorXd &coeffs) const {
        Eigen::VectorXd new_mc = reduce_zeros(coeffs);
        std::size_t Norder = new_mc.size() - 1;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Norder, Norder);
        // c_wrap wraps the first 0...Norder elements of the coefficient vector
        Eigen::Map<const Eigen::VectorXd> c_wrap(&(new_mc[0]), Norder);
        // First row
        A(0, 1) = 1;
        // Last row
        A.row(Norder - 1) = -c_wrap / (2.0*new_mc(Norder));
        A(Norder - 1, Norder - 2) += 0.5;
        // All the other rows
        for (int j = 1; j < Norder - 1; ++j) {
            A(j, j - 1) = 0.5;
            A(j, j + 1) = 0.5;
        }
        return A;
    }
    std::vector<double> ChebyshevExpansion::real_roots2(bool only_in_domain) const {
        //vector of roots to be returned
        std::vector<double> roots;

        auto N = m_c.size()-1;
        auto Ndegree_scaled = N*2;
        Eigen::VectorXd xscaled = get_CLnodes(Ndegree_scaled), yy = y_Clenshaw_xscaled(xscaled);
        
        // a,b,c can also be obtained by solving the matrix system:
        // [x_k^2, x_k, 1] = [b_k] for k in 1,2,3
        for (auto i = 0; i+2 < Ndegree_scaled+1; i += 2){
            const double &x_1 = xscaled[i + 0], &y_1 = yy[i + 0],
                         &x_2 = xscaled[i + 1], &y_2 = yy[i + 1],
                         &x_3 = xscaled[i + 2], &y_3 = yy[i + 2];
            double d = (x_3 - x_2)*(x_2 - x_1)*(x_3 - x_1);
            double a = ((x_3 - x_2)*y_1 - (x_3 - x_1)*y_2 + (x_2 - x_1)*y_3) / d;
            double b = (-(POW2(x_3) - POW2(x_2))*y_1 + (POW2(x_3) - POW2(x_1))*y_2 - (POW2(x_2) - POW2(x_1))*y_3) / d;
            double c = ((x_3 - x_2)*x_2*x_3*y_1 - (x_3 - x_1)*x_1*x_3*y_2 + (x_2 - x_1)*x_2*x_1*y_3) / d;

            // Discriminant of quadratic
            double D = b*b - 4*a*c;
            if (D >= 0) {
                double root1, root2;
                if (a == 0){
                    // Linear
                    root1 = -c/b;
                    root2 = -1000; // Something outside the domain; we are in scaled coordinates so this will definitely get rejected
                }
                else if (D == 0) { // Unlikely due to numerical precision
                    // Two equal real roots
                    root1 = -b/(2*a);
                    root2 = -1000; // Something outside the domain; we are in scaled coordinates so this will definitely get rejected
                }
                else {
                    // Numerically stable method for solving quadratic
                    // From https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
                    double sqrtD = sqrt(D);
                    if (b >= 0){
                        root1 = (-b - sqrtD)/(2*a);
                        root2 = 2*c/(-b-sqrtD);
                    }
                    else {
                        root1 = 2*c/(-b+sqrtD);
                        root2 = (-b+sqrtD)/(2*a);
                    }
                }
                bool in1 = inbetween(x_1, x_3, root1), in2 = inbetween(x_1, x_3, root2);
                const ChebyshevExpansion &e = *this;
                auto secant = [e](double a, double ya, double b, double yb, double yeps = 1e-14, double xeps = 1e-14) {
                    auto c = b - yb*(b - a) / (yb - ya);
                    auto yc = e.y_Clenshaw_xscaled(c);
                    for (auto i = 0; i < 50; ++i){
                        if (yc*ya > 0) {
                            a=c; ya=yc;
                        }
                        else {
                            b=c; yb=yc;
                        }
                        if (std::abs(b - a) < xeps) { break; }
                        if (std::abs(yc) < yeps){ break; }
                        c = b - yb*(b - a) / (yb - ya);
                        yc = e.y_Clenshaw_xscaled(c);
                    }
                    return c;
                };
                int Nroots_inside = static_cast<int>(in1) + static_cast<int>(in2);
                if (Nroots_inside == 2) {
                    // Split the domain at the midline of the quadratic, polish each root against the underlying expansion
                    double x_m = -b/(2*a), y_m = e.y_Clenshaw_xscaled(x_m);
                    root1 = secant(x_1, y_1, x_m, y_m);
                    root2 = secant(x_m, y_m, x_3, y_3);
                    // Rescale back into real-world values in [xmin,xmax] from [-1,1]
                    roots.push_back(((m_xmax - m_xmin)*root1 + (m_xmax + m_xmin)) / 2.0); 
                    roots.push_back(((m_xmax - m_xmin)*root2 + (m_xmax + m_xmin)) / 2.0);
                }
                else if(Nroots_inside == 1) {
                    root1 = secant(x_1, y_1, x_3, y_3);
                    roots.push_back(((m_xmax - m_xmin)*root1 + (m_xmax + m_xmin)) / 2.0);
                }
                else {}
            }
        }
        return roots;
    }
    std::vector<double> ChebyshevExpansion::real_roots(bool only_in_domain) const {
      //vector of roots to be returned
        std::vector<double> roots;
        Eigen::VectorXd new_mc = reduce_zeros(m_c);
        //if the Chebyshev polynomial is just a constant, then there are no roots
        //if a_0=0 then there are infinite roots, but for our purposes, infinite roots doesnt make sense
        if (new_mc.size()<=1){ //we choose <=1 to account for the case of no coefficients
          return roots; //roots is empty
        }

        //if the Chebyshev polynomial is linear, then the only possible root is -a_0/a_1
        //we have this else if block because eigen is not a fan of 1x1 matrices
        else if (new_mc.size()==2){
          double val_n11 = -new_mc(0)/new_mc(1);
          const bool is_in_domain = (val_n11 >= -1.0 && val_n11 <= 1.0);
          // Keep it if it is in domain, or if you just want all real roots
          if (!only_in_domain || is_in_domain) {
              // Rescale back into real-world values in [xmin,xmax] from [-1,1]
              double x = ((m_xmax - m_xmin)*val_n11 + (m_xmax + m_xmin)) / 2.0;
              roots.push_back(x);
          }
        }

        //this for all cases of higher order polynomials
        else{
          // The companion matrix is definitely lower Hessenberg, so we can skip the Hessenberg
          // decomposition, and get the real eigenvalues directly.  These eigenvalues are defined
          // in the domain [-1, 1], but it might also include values outside [-1, 1]
          Eigen::VectorXcd eigs = eigenvalues(companion_matrix(new_mc), /* balance = */ true);


          for (Eigen::Index i = 0; i < eigs.size(); ++i) {
              if (std::abs(eigs(i).imag() / eigs(i).real()) < 1e-15) {
                  double val_n11 = eigs(i).real();
                  const bool is_in_domain = (val_n11 >= -1.0 && val_n11 <= 1.0);
                  // Keep it if it is in domain, or if you just want all real roots
                  if (!only_in_domain || is_in_domain) {
                      // Rescale back into real-world values in [xmin,xmax] from [-1,1]
                      double x = ((m_xmax - m_xmin)*val_n11 + (m_xmax + m_xmin)) / 2.0;
                      roots.push_back(x);
                  }
              }
          }
        }
        return roots;

        //// The companion matrix is definitely lower Hessenberg, so we can skip the Hessenberg
        //// decomposition, and get the real eigenvalues directly.  These eigenvalues are defined
        //// in the domain [-1, 1], but it might also include values outside [-1, 1]
        //Eigen::VectorXd real_eigs = eigenvalues_upperHessenberg(companion_matrix().transpose(), /* balance = */ true);
        //
        //std::vector<double> roots;
        //for (Eigen::Index i = 0; i < real_eigs.size(); ++i){
        //    double val_n11 = real_eigs(i);
        //    const bool is_in_domain = (val_n11 >= -1.0 && val_n11 <= 1.0);
        //    // Keep it if it is in domain, or if you just want all real roots
        //    if (!only_in_domain || is_in_domain){
        //        // Rescale back into real-world values in [xmin,xmax] from [-1,1]
        //        double x = ((m_xmax - m_xmin)*val_n11 + (m_xmax + m_xmin)) / 2.0;
        //        roots.push_back(x);
        //    }
        //}
        //return roots;
    }
    double ChebyshevExpansion::monotonic_solvex(double y) {
        /*
        Function is known to be monotonic, so we can shortcut some of the solving steps used
        We don't use the eigenvalue method because it is too slow
        */
        auto& e = *this;
        auto secant = [e,y](double a, double ya, double b, double yb, double yeps = 1e-14, double xeps = 1e-14) {
            double c, yc;
            for (auto i = 0; i < 50; ++i) {
                c = b - yb * (b - a) / (yb - ya);
                yc = e.y_Clenshaw_xscaled(c)-y;
                if (yc * ya > 0) {
                    a = c; ya = yc;
                }
                else {
                    b = c; yb = yc;
                }
                if (std::abs(b - a) < xeps) { break; }
                if (std::abs(yc) < yeps) { break; }
            }
            return c;
        };
        // Determine if the function is monotonically increasing or decreasing
        auto ynodes = get_node_function_values();
        auto nodes = get_nodes_n11();
        bool increasing = ynodes[ynodes.size() - 1] > ynodes[0]; // Nodes go from 1 to -1 (stuck with this), but increasing says whether the value at the last *index* (x=-1) is greater than that of the first index (x=1).
        if (increasing) {
            if (y > ynodes[ynodes.size() - 1]) {
                throw std::invalid_argument("Argument is outside the range of the expansion");
            }
            if (y < ynodes[0]) {
                throw std::invalid_argument("Argument is outside the range of the expansion");
            }
        }
        else {
            if (y < ynodes[ynodes.size() - 1]) {
                throw std::invalid_argument("Argument is outside the range of the expansion");
            }
            if (y > ynodes[0]) {
                throw std::invalid_argument("Argument is outside the range of the expansion");
            }
        }

        // Interval bisection to find the Chebyshev-Lobatto nodes that bound the solution by bisection
        int N = static_cast<int>(ynodes.size());
        Eigen::Index i = (increasing) ? get_increasingleftofval(ynodes, y, N) : get_decreasingleftofval(ynodes, y, N);
        auto xscaled = secant(nodes(i), ynodes(i)-y, nodes(i + 1), ynodes(i + 1)-y);
        return unscale_x(xscaled);
    }

    std::vector<ChebyshevExpansion> ChebyshevExpansion::subdivide(std::size_t Nintervals, const std::size_t Norder) const {

        if (Nintervals == 1) {
            return std::vector<ChebyshevExpansion>(1, *this);
        }

        std::vector<ChebyshevExpansion> segments;
        double deltax = (m_xmax - m_xmin) / (Nintervals - 1);

        // Vector of values in the range [-1,1] as roots of a high-order Chebyshev
        double NN = static_cast<double>(Norder);
        Eigen::VectorXd xpts_n11 = (Eigen::VectorXd::LinSpaced(Norder + 1, 0, NN)*EIGEN_PI/NN).array().cos();

        for (std::size_t i = 0; i < Nintervals - 1; ++i) {
            double xmin = m_xmin + i*deltax, xmax = m_xmin + (i + 1)*deltax;
            Eigen::VectorXd xrealworld = ((xmax - xmin)*xpts_n11.array() + (xmax + xmin)) / 2.0;
            segments.push_back(factoryf(Norder, y(xrealworld), xmin, xmax));
        }
        return segments;
    }
    std::vector<double> ChebyshevExpansion::real_roots_intervals(const std::vector<ChebyshevExpansion> &segments, bool only_in_domain) {
        std::vector<double> roots;
        for (auto &seg : segments) {
            const auto segroots = seg.real_roots(only_in_domain);
            roots.insert(roots.end(), segroots.cbegin(), segroots.cend());
        }
        return roots;
    }
    std::vector<double> ChebyshevExpansion::real_roots_approx(long Npoints)
    {
        std::vector<double> roots;
        // Vector of values in the range [-1,1] as roots of a high-order Chebyshev
        Eigen::VectorXd xpts_n11 = (Eigen::VectorXd::LinSpaced(Npoints + 1, 0, Npoints)*EIGEN_PI / Npoints).array().cos();
        // Scale values into real-world values
        Eigen::VectorXd ypts = y_recurrence_xscaled(xpts_n11);
        // Eigen::MatrixXd buf(Npoints+1, 2); buf.col(0) = xpts; buf.col(1) = ypts; std::cout << buf << std::endl;
        for (size_t i = 0; i < Npoints - 1; ++i) {
            // The change of sign guarantees at least one root between indices i and i+1
            double y1 = ypts(i), y2 = ypts(i + 1);
            bool signchange = (std::signbit(y1) != std::signbit(y2));
            if (signchange) {
                double xscaled = xpts_n11(i);

                // Fit a quadratic given three points; i and i+1 bracket the root, so need one more constraint
                // i0 is the leftmost of the three indices that will be used; when i == 0, use
                // indices i,i+1,i+2, otherwise i-1,i,i+1
                size_t i0 = (i >= 1) ? i - 1 : i;
                Eigen::Vector3d r;
                r << ypts(i0), ypts(i0 + 1), ypts(i0 + 2);
                Eigen::Matrix3d A;
                for (std::size_t irow = 0; irow < 3; ++irow) {
                    double _x = xpts_n11(i0 + irow);
                    A.row(irow) << _x*_x, _x, 1;
                }
                // abc holds the coefficients a,b,c for y = a*x^2 + b*x + c
                Eigen::VectorXd abc = A.colPivHouseholderQr().solve(r);
                double a = abc[0], b = abc[1], c = abc[2];

                // Solve the quadratic and find the root you want
                double x1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
                double x2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
                bool x1_in_range = is_in_closed_range(xpts_n11(i), xpts_n11(i + 1), x1);
                bool x2_in_range = is_in_closed_range(xpts_n11(i), xpts_n11(i + 1), x2);

                // Double check that only one root is within the range
                if (x1_in_range && !x2_in_range) {
                    xscaled = x1;
                }
                else if (x2_in_range && !x1_in_range) {
                    xscaled = x2;
                }
                else {
                    xscaled = 1e99;
                }

                // Rescale back into real-world values
                double x = ((m_xmax - m_xmin)*xscaled + (m_xmax + m_xmin)) / 2.0;
                roots.push_back(x);
            }
            else {
                // TODO: locate other roots based on derivative considerations
            }
        }
        return roots;
    }

    /// Chebyshev-Lobatto nodes \f$ \cos(\pi j/N), j = 0,..., N \f$ in the range [-1,1]
    Eigen::VectorXd ChebyshevExpansion::get_nodes_n11() {
        std::size_t N = m_c.size()-1;
        return CLnodes_library.get(N);
    }
    /// Chebyshev-Lobatto nodes \f$\cos(\pi j/N), j = 0,..., N \f$ mapped to the range [xmin, xmax]
    Eigen::VectorXd ChebyshevExpansion::get_nodes_realworld() {
        return ((m_xmax - m_xmin)*get_nodes_n11().array() + (m_xmax + m_xmin))*0.5;
    }
    /// Values of the function at the Chebyshev-Lobatto nodes 
    Eigen::VectorXd ChebyshevExpansion::get_node_function_values() const{
        if (m_nodal_value_cache.size() > 0) {
            return m_nodal_value_cache;
        }
        else {
            std::size_t N = m_c.size() - 1;
            return u_matrix_library.get(N) * m_c;
        }
    }
    ChebyshevExpansion ChebyshevExpansion::factoryf(const std::size_t N, const Eigen::VectorXd &f, const double xmin, const double xmax) {
        // Step 3: Get coefficients for the L matrix from the library of coefficients
        const Eigen::MatrixXd &L = l_matrix_library.get(N);
        // Step 4: Obtain coefficients from vector - matrix product
        return ChebyshevExpansion(L*f, xmin, xmax);
    }
    ChebyshevExpansion ChebyshevExpansion::factoryfFFT(const std::size_t N, const Eigen::VectorXd& f, const double xmin, const double xmax) {

        Eigen::VectorXd valsUnitDisc(2 * f.size() - 2);
        // Starting at x = 1, going to -1, then the same nodes, not including x=-1 and x=1, in the opposite order
        valsUnitDisc.head(f.size()) = f;
        valsUnitDisc.tail(f.size() - 2) = f.reverse().segment(1, f.size() - 2);
        
        Eigen::FFT<double> fft;
        Eigen::VectorXcd FourierCoeffs(2 * f.size() - 2);
        fft.fwd(FourierCoeffs, valsUnitDisc);
        auto n = f.size() - 1;
        Eigen::ArrayXd ChebCoeffs = FourierCoeffs.real().head(n+1)/n;
        ChebCoeffs[0] /= 2;
        ChebCoeffs[ChebCoeffs.size()-1] /= 2;

        return ChebyshevExpansion(ChebCoeffs, xmin, xmax);
    }
    ChebyshevExpansion ChebyshevExpansion::from_powxn(const std::size_t n, const double xmin, const double xmax) {
        if (xmin != -1) {
            throw std::invalid_argument("xmin must be -1");
        }
        if (xmax != 1) {
            throw std::invalid_argument("xmax must be 1");
        }
        Eigen::VectorXd c = Eigen::VectorXd::Zero(n + 1);
        for (std::size_t k = 0; k <= n / 2; ++k) {
            std::size_t index = n - 2 * k;
            double coeff = binomialCoefficient(static_cast<double>(n), static_cast<double>(k));
            if (index == 0) {
                coeff /= 2.0;
            }
            c(index) = coeff;
        }
        return pow(2, 1-static_cast<int>(n))*ChebyshevExpansion(c, xmin, xmax);
    }
    ChebyshevExpansion ChebyshevExpansion::deriv(std::size_t Nderiv) const {
        // See Mason and Handscomb, p. 34, Eq. 2.52
        // and example in https ://github.com/numpy/numpy/blob/master/numpy/polynomial/chebyshev.py#L868-L964
        vectype c = m_c;
        for (std::size_t deriv_counter = 0; deriv_counter < Nderiv; ++deriv_counter) {
            std::size_t N = c.size() - 1, ///< Order of the expansion
                        Nd = N - 1; ///< Order of the derivative expansion
            vectype cd(N);
            for (std::size_t r = 0; r <= Nd; ++r) {
                cd(r) = 0;
                for (std::size_t k = r + 1; k <= N; ++k) {
                    // Terms where k-r is odd have values, otherwise, they are zero
                    if ((k - r) % 2 == 1) {
                        cd(r) += 2*k*c(k);
                    }
                }
                // The first term with r = 0 is divided by 2 (the single prime in Mason and Handscomb, p. 34, Eq. 2.52)
                if (r == 0) {
                    cd(r) /= 2;
                }
                // Rescale the values if the range is not [-1,1].  Arrives from the derivative of d(xreal)/d(x_{-1,1})
                cd(r) /= (m_xmax-m_xmin)/2.0;
            }
            if (Nderiv == 1) {
                return ChebyshevExpansion(std::move(cd), m_xmin, m_xmax);
            }
            else{
                c = cd;
            }
        }
        return ChebyshevExpansion(std::move(c), m_xmin, m_xmax);
    };
    ChebyshevExpansion ChebyshevExpansion::integrate(std::size_t Nintegral) const {
        // See Mason and Handscomb, p. 33, Eq. 2.44 & 2.45
        // and example in https ://github.com/numpy/numpy/blob/master/numpy/polynomial/chebyshev.py#L868-L964
        if (Nintegral != 1) { throw std::invalid_argument("Only support one integral for now"); }
        vectype c(m_c.size() + 1);
        double width = m_xmax - m_xmin;
        for (auto i = 1; i < m_c.size()+1; ++i) {
            if (i == 1) {
                // This special case is needed because the prime on the summation in Mason indicates the first coefficient 
                // is to be divided by two
                c[i] = (2*m_c[i - 1] - m_c[i + 1]) / (2 * i);
            }
            else if (i + 1 > m_c.size()-1) {
                c[i] = (m_c[i - 1]) / (2 * i);
            }
            else {
                c[i] = (m_c[i - 1] - m_c[i + 1]) / (2 * i);
            }
        }
        c(0) = 0; // This is the arbitrary constant;
        c *= width / 2;
        return ChebyshevExpansion(std::move(c), m_xmin, m_xmax);
    }

    Eigen::VectorXd eigenvalues_upperHessenberg(const Eigen::MatrixXd &A, bool balance){
        Eigen::VectorXd roots(A.cols());
        Eigen::RealSchur<Eigen::MatrixXd> schur;

        if (balance) {
            Eigen::MatrixXd Abalanced, D;
            balance_matrix(A, Abalanced, D);
            schur.computeFromHessenberg(Abalanced, Eigen::MatrixXd::Zero(Abalanced.rows(), Abalanced.cols()), false);
        }
        else {
            schur.computeFromHessenberg(A, Eigen::MatrixXd::Zero(A.rows(), A.cols()), false);
        }

        const Eigen::MatrixXd &T = schur.matrixT();
        Eigen::Index j = 0;
        for (int i = 0; i < T.cols(); ++i) {
            if (i+1 < T.cols()-1 && std::abs(T(i+1,i)) > DBL_EPSILON){
                // Nope, this is a 2x2 block, keep moving
                i += 1;
            }
            else{
                // This is a 1x1 block, keep this (real) eigenvalue
                roots(j) = T(i, i);
                j++;
            }
        }
        roots.conservativeResize(j-1);
        return roots;
    }

    Eigen::VectorXcd eigenvalues(const Eigen::MatrixXd &A, bool balance) {
        if (balance) {
            Eigen::MatrixXd Abalanced, D;
            balance_matrix(A, Abalanced, D);
            return Abalanced.eigenvalues();
        }
        else {
            return A.eigenvalues();
        }
    }

}; /* namespace ChebTools */
