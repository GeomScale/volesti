/**
 * @file BotevAlgorithm.cpp
 *
 * @brief A collection of functions to deal with the truncated
 * univariate and multivariate normal distributions.
 *
 * @author Zdravko Botev (original R version)
 * @author Shinsuke Mori (Hiroshima University, port to C++)
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Mutsuo Saito (little bit modification for publish)
 *
 * Copyright (C) 2015 Zdravko I. Botev.
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 *
 * @note ported from
 * https://cran.r-project.org/package=TruncatedNormal
 * under License of GPL-3(https://cran.r-project.org/web/licenses/GPL-3)
 */

#ifndef BETOVENALGORITHMTWO_H
#define BETOVENALGORITHMTWO_H

#include "pnorm2.h" // for log_pnorm(), pnorm() and pnorm5()
#include "kahan.hpp" // for cholperm() and so on
#include "BotevAlgorithm.h"
#include "DigitalNet2.h" // for mvnprqmc()

#include <vector> // for cholperm()
#include <math.h> // for cholperm and so on
#include <random> // for random_device
#include <iostream>
#include <iomanip>
#include <random>
#include <boost/numeric/ublas/lu.hpp>
#include <stdio.h>

//#define DEBUG 1

using namespace std;

typedef boost::numeric::ublas::vector<double> Vectord;
typedef boost::numeric::ublas::matrix<double> Matrixd;

/**
 * Unnamed Namespace for File Scope.
 */
namespace {
    const double EPS = 1e-10;
    const double one_over_sqrt2pi
    = 0.398942280401432677939946059934381868475858631164934657665;
    //const int PROB = 99; // 95, 99, 999, 9999 is available

#define MVN_MAGIC UINT64_C(0x2d221e655a0f6651)

    struct mvn_parameter_name {
        std::string name;
        std::string abb;
        std::string construction;
    };

    const mvn_parameter_name mvn_parameter_name_data[] = {
        {"Miwa", "miwa", "l = -inf, u = 0, matSigma = {2s/(s+1), s/(s+1)}"},
        {"Random1", "rnd1", "l = -inf, u = 0, matSigma: random"},
        {"Random2", "rnd2", "l, u, matSigma: random"},
        {"Random3", "rnd3", "l, u, matSigma: random"}
    };

    const uint32_t mvn_parameter_name_data_size = 4;

    struct mvn_file_header {
        uint32_t s;
        uint32_t pos;
    };

    class MVNFileData {
    public:
        MVNFileData(uint32_t s) {
            this->s = s;
            this->lower = new double[s];
            this->upper = new double[s];
            this->sigma = new double[s * s];
        }
        ~MVNFileData() {
            delete[] lower;
            delete[] upper;
            delete[] sigma;
        }
        uint32_t s;
        double * lower;
        double * upper;
        double * sigma;
    };

    /**
     * Cholesky Decomposition
     *
     * If sigma is not pisitive definite, this function returns false,
     * otherwise, returns true and set L.
     *
     * Note: L should have the same size as sigma.
     * sigma -> L
     * @param[in] s number of rows and columns.
     * @param[in] sigma s * s symmetric real matrix.
     * @param[out] L the result of Cholesky Decomposition of sigma.
     * @return true if sigma is positive definite
     */
    bool choleskyDecomposition(int s,
                               const boost::numeric::ublas::matrix<double>& sigma,
                               boost::numeric::ublas::matrix<double>& L)
    {
        for (int i = 0; i < s; ++i) {
            for (int j = 0; j < s; ++j) {
                if ( fabs(sigma(i, j) - sigma(j, i)) > 1e-10 ) {
                    //cout << "Error: sigma is not symmetric." << endl;
                    //throw runtime_error("sigma is not symmetric.");
                    return false;
                }
                if (j < i) {
                    Kahan sum;
                    for (int k = 0; k < j; ++k) {
                        sum.add(L(i, k) * L(j, k));
                    }
                    L(i, j) = (sigma(i, j) - sum.get()) / L(j, j);
                } else if (j == i) {
                    Kahan sum;
                    for (int k = 0; k < j; ++k) {
                        sum.add(L(i, k) * L(i, k));
                    }
                    double tmp = sigma(i, j) - sum.get();
                    if (tmp <= 0.0) {
                        //cout << "Error: sigma is not positive-definite."
                        //<< endl;
                        //throw runtime_error(
                        //"sigma is not positive-definite.");
                        return false;
                    }
                    L(i, j) = sqrt(tmp);
                } else {
                    L(i, j) = 0.0;
                }
            }
        }
        return true;
    }

    /*
      Functions used in Botev Algorithm
    */
    template <typename T> void swap(T *ptr1, T *ptr2) {
        T tmp = *ptr1;
        *ptr1 = *ptr2;
        *ptr2 = tmp;
    }

    /**
     * computes ln(P(a<Z<b))
     * where Z~N(0,1) very accurately for any 'a', 'b'
     */
    inline double lnNpr(double lower, double upper)
    {
        if (0.0 < lower) {
            double pl = log_pnorm(-lower);
            double pu = log_pnorm(-upper);
            return pl + log1p(-exp(pu - pl));
        } else if (upper < 0.0) {
            double pl = log_pnorm(lower);
            double pu = log_pnorm(upper);
            return pu + log1p(-exp(pl - pu));
        } else {
            double pl = pnorm(lower);
            double pu = pnorm(-upper);
            return log1p(- pl - pu);
        }
    }

    /**
     *
     *
     */
    inline int argminN01(std::vector<double> lower,
                         std::vector<double> upper)
    {
        std::vector<double> v;
        const int size = lower.size();
        for (int i = 0; i < size; ++i) {
            v.push_back(lnNpr(lower[i], upper[i]));
        }
        return distance(v.begin(), min_element(v.begin(), v.end()));
    }

    /**
     * define Q function
     */
    inline double qfun(double x) {
        return exp(0.5 * x * x + pnorm5(x, 0.0, 1.0, 0, 1));
    }

    /**
     * computes with precision the quantile function
     * of the standard normal distribution,
     * truncated to the interval [l,u];
     * normq assumes 0<l<u and 0<p<1
     * lower:l
     * upper:u
     * @param p
     * @param lower
     * @param upper
     * @return
     */
    inline double normq(double p, double lower, double upper)
    {
        using namespace std;
        if (lower > 100000.0) {
            return sqrt(lower * lower - 2.0
                        * log(1.0 + p
                              * expm1(lower * lower * 0.5
                                      - upper * upper * 0.5)));
        }
        // newton(p, lower, upper);
        const double ql = qfun(lower);
        double qu = 0.0;
        if ( !isinf(upper) )
        {
            qu = qfun(upper);
        }
        lower *= lower;
        upper *= upper;
        double x = sqrt(lower - 2.0
                        * log(1.0 + p * expm1(lower * 0.5 - upper * 0.5)));
        double err = INFINITY;
        while ( err > EPS ) {
            double del = -qfun(x) + (1.0 - p)
                * exp(0.5 * (x * x - lower))*ql + p
                * exp(0.5 * (x * x - upper)) * qu;
            x -= del;
            err = fabs(del);
        }
        return x;
    }

    /*
     *   Here is original code: https://gist.github.com/kmpm/1211922/
     * Original C++ implementation found at
     * http://www.wilmott.com/messageview.cfm?catid=10&threadid=38771
     * C# implementation found at
     * http://weblogs.asp.net/esanchez/archive/2010/07/29/a-quick-and-dirty-implementation-of-excel-norminv-function-in-c.aspx
     */
    inline double Phiinv(double p)
    {
        if ( p == 0.0 ) {
            return -INFINITY;
        }
        if ( p == 1.0 ) {
            return INFINITY;
        }
        double q = p - 0.5;
        double ret;
        if (fabs(q) <= 0.425) {
            double r = 0.180625 - q * q;
            ret = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105)
                            * r + 67265.770927008700853)
                           * r + 45921.953931549871457)
                          * r + 13731.693765509461125)
                         * r + 1971.5909503065514427)
                        * r + 133.14166789178437745)
                       * r + 3.387132872796366608)
                / (((((((r * 5226.495278852854561 + 28729.085735721942674)
                        * r + 39307.89580009271061)
                       * r + 21213.794301586595867)
                      * r + 5394.1960214247511077)
                     * r + 687.1870074920579083)
                    * r + 42.313330701600911252)
                   * r + 1);
        } else {
            double r;
            if (q > 0.0) {
                r = 1.0 - p;
            } else {
                r = p;
            }
            r = sqrt(-log(r));
            if ( r <= 5 ) {
                r -= 1.6;
                ret = (((((((r * 7.7454501427834140764e-4
                             + .0227238449892691845833)
                            * r + .24178072517745061177)
                           * r + 1.27045825245236838258)
                          * r + 3.64784832476320460504)
                         * r + 5.7694972214606914055)
                        * r + 4.6303378461565452959)
                       * r + 1.42343711074968357734)
                    / (((((((r *1.05075007164441684324e-9
                             + 5.475938084995344946e-4)
                            * r + .0151986665636164571966)
                           * r + .14810397642748007459)
                          * r + .68976733498510000455)
                         * r + 1.6763848301838038494)
                        * r + 2.05319162663775882187)
                       * r + 1);
            } else {
                r -= 5.0;
                ret = (((((((r * 2.01033439929228813265e-7
                             + 2.71155556874348757815e-5)
                            * r + .0012426609473880784386)
                           * r + .026532189526576123093)
                          * r + .29656057182850489123)
                         * r + 1.7848265399172913358)
                        * r + 5.4637849111641143699)
                       * r + 6.6579046435011037772)
                    / (((((((r * 2.04426310338993978564e-15
                             + 1.4215117583164458887e-7)
                            * r + 1.8463183175100546818e-5)
                           * r + 7.868691311456132591e-4)
                          * r + .0148753612908506148525)
                         * r + .13692988092273580531)
                        * r + .59983220655588793769)
                       * r + 1);
            }
            if ( q < 0.0 ) {
                ret = -ret;
            }
        }
        return ret;
    }

    /*
     * computes with precision the quantile function
     * of the standard normal distribution,
     * truncated to the interval [l,u], using erfcinv.
     * Phiinv.R in TruncatedNormal_1.0.tar.gz
     */
    inline double Phinv(double p, double lower, double upper)
    {
        if (upper < 0.0) {
            return -Phiinv(pnorm(-lower)+(pnorm(-upper) - pnorm(-lower)) * p);
        }
        return Phiinv(pnorm(lower)+(pnorm(upper) - pnorm(lower)) * p);
    }

    /**
     * normal quantile function with precision
     * computes with tail-precision the quantile function
     * of the standard normal distribution at 0<=p<=1,
     * and truncated to the interval [l,u];
     * Inf values for vectors 'l' and 'u' accepted;
     *
     * * Example 1:
     * # Suppose you wish to simulate a random variable
     * # 'Z' from the non-standard Gaussian N(m,s^2)
     * # conditional on l<Z<u. First compute
     *  m=1;l=10;u=20;s=1;
     *  X=norminvp(runif(1),(l-m)/s,(u-m)/s); # and then set
     *  Z=m+s*X
     *
     * Reference:
     * Z. I. Botev (2015),
     * "The Normal Law Under Linear Restrictions:
     *  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
     *
     * @param p
     * @param lower
     * @param upper
     * @return
     */
    inline double norminvp(double p, double lower, double upper)
    {
        if (p == 0.0) {
            return lower;
        }
        if (p == 1.0) {
            return upper;
        }
        const double TRESHOLD = 35.0;
        if (TRESHOLD < lower) {
            return normq(p, lower, upper);
        }
        if (upper < -TRESHOLD) {
            return -normq(1.0 - p, lower, upper);
        }
        return Phinv(p, lower, upper);
    }

    /**
     * rescale lower, upper and L
     *
     * @param[in] s dimension
     * @param[in/out] lower lower bound
     * @param[in/out] upper upper bound
     * @param[in/out] L Cholesky Decomposition of sigma
     */
    void rescale(int s,
                 Vectord& lower,
                 Vectord& upper,
                 Matrixd& L)
    {
        for (int i = 0; i < s; ++i) {
            lower(i) /= L(i, i);
            upper(i) /= L(i, i);
            for (int j = 0; j < i; ++j) {
                L(i, j) /= L(i, i);
            }
            L(i, i) = 0.0;
        }
    }

    /**
     * implements grad_psi(x) to find optimal exponential twisting;
     * assume scaled 'L' with zero diagonal;
     *
     * @param[in] s dimension
     * @param[in] y
     * @param[in/out] L
     * @param[in] lower lower bound:l
     * @param[in] upper upper bound:u
     * @param[out] Jacobian
     * @param[out] grad the result
     */
    void gradpsi(int s,
                 const Vectord& y,
                 const Matrixd& L,
                 const Vectord& lower,
                 const Vectord& upper,
                 Matrixd& Jacobian,
                 Vectord& grad)
    {
        Vectord c(s);
        c.clear(); // clear() set all element zero
        // 何やってんの？
        //Vectord x = c;
        //Vectord mu = c;
        Vectord x(s);
        Vectord mu(s);
        // 本当か?
        x(s - 1) = c(s - 1);
        mu(s - 1) = c(s - 1);
        for (int i = 0; i < s-1; ++i) {
            x(i) = y(i);
            mu(i) = y(s-1+i);
        }
        c(0) = 0.0;
        for (int i = 1; i < s; ++i) {
            Kahan sum;
            for (int k = 0; k < i; ++k) {
                sum.add(L(i, k) * x(k));
            }
            c(i) = sum.get();
        }
        Vectord lt = lower - mu - c;
        Vectord ut = upper - mu - c;
        Vectord w(s);
        for (int i = 0; i < s; ++i) {
            w(i) = lnNpr(lt(i), ut(i));
        }
        Vectord P(s);
        Vectord pl(s);
        Vectord pu(s);
        for (int i = 0; i < s; ++i) {
            pl(i) = exp(-0.5 * lt(i) * lt(i) - w(i)) * one_over_sqrt2pi;
            pu(i) = exp(-0.5 * ut(i) * ut(i) - w(i)) * one_over_sqrt2pi;
            P(i) = pl(i) - pu(i);
        }
        Vectord dfdx = mu - trans(prod(trans(P), L));
        dfdx.resize(s-1);
        Vectord dfdm = - mu + x - P;
        dfdm.resize(s-1);
        grad.resize(2*s-2);
        for (int i = 0; i < s-1; ++i) {
            grad(i) = dfdx(i);
            grad(s-1+i) = dfdm(i);
        }
        Vectord dP(s);
        for (int i = 0; i < s; ++i) {
            if ( isinf(lt(i)) ) {
                lt(i) = 0.0;
            }
            if ( isinf(ut(i)) ) {
                ut(i) = 0.0;
            }
            dP(i) = - P(i) * P(i) + lt(i) * pl(i) - ut(i) * pu(i);
        }
        Matrixd DL(s, s);
        for (int i = 0; i < s; ++i) {
            for (int j = 0; j < i; ++j) {
                DL(i, j) = dP(i) * L(i, j);
            }
            for (int j = i; j < s; ++j) {
                DL(i, j) = 0.0;
            }
        }
        Matrixd mx = DL;
        for (int i = 0; i < s; ++i) {
            mx(i, i) -= 1.0;
        }
        mx.resize(s-1, s-1);
        Matrixd xx = prod(trans(L), DL);
        xx.resize(s-1, s-1);
        Jacobian.resize(2*s-2, 2*s-2);
        Jacobian.clear();
        for (int j = 0; j < s-1; ++j) {
            for (int i = 0; i < s-1; ++i) {
                Jacobian(i, j) = xx(i, j);
            }
            for (int i = s-1; i < 2*s-2; ++i) {
                Jacobian(i, j) = mx(i-(s-1), j);
            }
        }
        for (int j = s-1; j < 2*s-2; ++j) {
            for (int i = 0; i < s-1; ++i) {
                Jacobian(i, j) = Jacobian(j, i);
            }
        }
        for (int i = s-1; i < 2*s-2; ++i) {
            for (int j = s-1; j < 2*s-2; ++j) {
                if ( i == j ) {
                    Jacobian(i, j) = 1.0 + dP(i-(s-1));
                } else {
                    Jacobian(i, j) = 0.0;
                }
            }
        }
    }

    /**
     *  Computes permuted lower Cholesky factor L for Sig
     *  by permuting integration limit vectors l and u.
     *  Outputs perm, such that Sig(perm,perm)=L%*%t(L).
     *
     * Reference:
     *  Gibson G. J., Glasbey C. A., Elston D. A. (1994),
     *  "Monte Carlo evaluation of multivariate normal integrals and
     *  sensitivity to variate ordering",
     *  In: Advances in Numerical Methods and Applications, pages 120--126
     *
     * cholperm <-
     *  function( Sig, l, u )
     *
     * @param[in] s dimension
     * @param[in/out] sigma variant/covariant matrix:Sig
     * @param[in/out] lower lower integration limit:l
     * @param[in/out] upper upper integration limit:u
     * @param[out] perm permutation information array of size s.
     * @param[out] L permuted Cholesky factor for sigma
     */
    void cholPerm(uint32_t s,
                  Matrixd& sigma,
                  Vectord& lower,
                  Vectord& upper,
                  int perm[],
                  Matrixd& L)
    {
#if defined(DEBUG)
        cout << "in cholPerm" << endl;
#endif
        //int s = mvn.s;
        for (uint32_t i = 0; i < s; ++i) {
            perm[i] = i;
        }
        L.resize(s, s);
        L.clear();
        double z[s];
        for (uint32_t i = 0; i < s; ++i) {
            z[i] = 0.0;
        }
        std::vector<double> tl;
        std::vector<double> tu;
        for (uint32_t j = 0; j < s; ++j) {
            tl.clear();
            tu.clear();
            for (uint32_t i = j; i < s; ++i) {
                Kahan sum_stmp;
                Kahan sum_cols;
                for (uint32_t k = 0; k < j; ++k) {
                    sum_stmp.add(L(i, k) * L(i, k));
                    sum_cols.add(L(i, k) * z[k]);
                }
                double stmp = sqrt(sigma(i, i) - sum_stmp.get());
                double cols = sum_cols.get();
                tl.push_back((lower(i) - cols)/stmp);
                tu.push_back((upper(i) - cols)/stmp);
            }
#if defined(DEBUG)
        cout << "before argminN01" << endl;
#endif
            int k = argminN01(tl, tu) + j;
            //flip(j, k);
            // start flip
            //int s = mvn.s;
            // l, u, perm のj番目とk番目の入れ替え
#if defined(DEBUG)
        cout << "after argminN01" << endl;
#endif
            swap(&lower(j), &lower(k));
            swap(&upper(j), &upper(k));
            swap(&perm[j], &perm[k]);
            for (uint32_t i = 0; i < s; ++i) {
                // Sigmaのj行とk行の入れ替え
                swap(&sigma(j, i), &sigma(k, i));
                // Lのj行とk行の入れ替え(列の入れ替えはない)
                swap(&L(j, i), &L(k, i));
            }
            for (uint32_t i = 0; i < s; ++i) {
                // Sigmaのj列とk列の入れ替え
                swap(&sigma(i, j), &sigma(i, k));
            }
            // end flip
            double stmp = sigma(j, j);
            Kahan sum;
            for (uint32_t k = 0; k < j; ++k) {
                sum.add(L(j, k) * L(j, k));
            }
            stmp -= sum.get();
            if ( stmp < 0.0 ) {
                cerr << "Error: Sigma is not positive semi-definite" << endl;
                throw runtime_error("sigma is not positive difinite.");
            }
            L(j, j) = sqrt(stmp);
            for (uint32_t i = j+1; i < s; ++i) {
                sum.clear();
                for (uint32_t k = 0; k < j; ++k) {
                    sum.add(L(i, k) * L(j, k));
                }
                L(i, j) = (sigma(i, j) - sum.get()) / L(j, j);
            }
            sum.clear();
            for (uint32_t k = 0; k <= j; ++k) {
                sum.add(L(j, k) * z[k]);
            }
            double tltmp = (lower(j) - sum.get()) / L(j, j);
            double tutmp = (upper(j) - sum.get()) / L(j, j);
            double wtmp = log(pnorm(tutmp) - pnorm(tltmp));
            z[j] = (exp(-0.5 * tltmp * tltmp -wtmp)
                    - exp(-0.5 * tutmp * tutmp -wtmp)) * one_over_sqrt2pi;
        }
    }

    double upbnd(uint32_t s,
                 const Vectord& lower,
                 const Vectord& upper,
                 const Matrixd& L,
                 const Vectord& xmu)
    {
        //int s = mvn.s;
        double x[s];
        double mu[s];
        for (uint32_t i = 0; i < s-1; ++i) {
            x[i] = xmu(i);
            mu[i] = xmu(s-1+i);
        }
        x[s-1] = 0.0;
        mu[s-1] = 0.0;
        Kahan ret;
        for (uint32_t i = 0; i < s; ++i) {
            Kahan sum;
            for (uint32_t j = 0; j < i; ++j) {
                sum.add(L(i, j) * x[j]);
            }
            ret.add(lnNpr(lower(i) - mu[i] - sum.get(), upper(i)
                          - mu[i] - sum.get())
                    + 0.5 * mu[i] * mu[i] - x[i] * mu[i]);
        }
        return exp(ret.get());
    }


    /*
     * non-linear equation solver
     * This method uses
     * boost::numeric::ublas::lu_factorize and
     * boost::numeric::ublas::lu_substitute.
     * Jacobian and grad should be boost::numeric::ublas::vector<double>
     */
    bool nleq(uint32_t s,
              const Vectord& lower,
              const Vectord& upper,
              const Matrixd& L,
              Vectord& xmu)
    {
#if defined(DEBUG)
        cout << "in nleq" << endl;
#endif
        //int s = mvn.s;
        double err = INFINITY;
        int iter = 0;
        Matrixd Jacobian;
        while ( err > EPS ) {
            //gradpsi(s, xmu, L, mvn.lower, mvn.upper, Jacobian);
            //Vectord del = grad;
            Vectord grad;
#if defined(DEBUG)
            cout << "before gradpsi" << endl;
#endif
            gradpsi(s, xmu, L, lower, upper, Jacobian, grad);
#if defined(DEBUG)
            cout << "after gradpsi" << endl;
#endif
            Vectord del = grad;

#if defined(DEBUG)
            cout << "before permutation_matrix" << endl;
#endif
            boost::numeric::ublas::permutation_matrix<> pm(2*s-2);
#if defined(DEBUG)
            cout << "before lu_factorize" << endl;
#endif
            boost::numeric::ublas::lu_factorize(Jacobian, pm);
#if defined(DEBUG)
            cout << "before lu_substitute" << endl;
            cout << "Jacobian:" << endl;
            for (uint32_t i = 0; i < Jacobian.size1(); i++) {
                for (uint32_t j = 0; j < Jacobian.size2(); j++) {
                    cout << Jacobian(i, j) << "\t";
                }
                cout << endl;
            }
            cout << endl;
            cout << "pm:" << endl;
            for (uint32_t i = 0; i < pm.size(); i++) {
                cout << pm(i) << "\t";
            }
            cout << endl;
            cout << "del:" << endl;
            for (uint32_t i = 0; i < del.size(); i++) {
                cout << del(i) << "\t";
            }
#endif
            boost::numeric::ublas::lu_substitute(Jacobian, pm, del);
#if defined(DEBUG)
            cout << "after lu_substitute" << endl;
#endif
            xmu += del;
            Kahan sum;
            for (uint32_t i = 0; i < 2*s-2; ++i) {
                sum.add(grad(i) * grad(i));
            }
            err = sum.get();
            iter++;
            if ( iter > 150 ) {
                cerr << "Covariance matrix is ill-conditioned and "
                     << "method failed." << endl;
                //throw runtime_error("Covariance matrix is ill-conditioned.");
                return false;
            }
        }
        return true;
    }


    int read_mvn(string& filename, MVNFileData& mvn)
    {
        uint32_t s = mvn.s;
        const char * mode = "rb";
        size_t count;
        FILE *fp = fopen(filename.c_str(), mode);
        //rewind(fp);
        uint64_t dmy;
        count = fread(&dmy, sizeof(uint64_t), 1, fp);
        if (count != 1) {
            cout << "fail to read magic number" << endl;
            return -1;
        }
        if (dmy != MVN_MAGIC) {
            cout << "magic number mismatch" << endl;
            return -1;
        }
        mvn_file_header header;
        for (;;) {
            count = fread(&header, sizeof(mvn_file_header), 1, fp);
            if (count != 1) {
                cout << "fail to read header s = "
                     << dec << s << endl;
                return -1;
            }
            if (header.s > s) {
                cout << "header s = " << dec << header.s << endl;
                return -1;
            }
            //cout << "header s = " << dec << header.s << endl;
            if (header.s == s) {
                break;
            }
        }
        fseek(fp, header.pos, SEEK_SET);
        count = fread(mvn.lower, sizeof(double), s, fp);
        if (count != s) {
            cout << "fail to read lower data" << endl;
            return -1;
        }
        count = fread(mvn.upper, sizeof(double), s, fp);
        if (count != s) {
            cout << "fail to read upper data" << endl;
            return -1;
        }
        size_t size = s * s;
        count = fread(mvn.sigma, sizeof(double), size, fp);
        if (count != size) {
            cout << "fail to read sigma data" << endl;
            return -1;
        }
        return 0;
    }

    int read_mvn(int id, int s, MVNFileData& mvn)
    {
        string filename = mvn_parameter_name_data[id].name;
        mvn.s = s;
        return read_mvn(filename, mvn);
    }
}

namespace MCQMCIntegration {
    class BotevAlgorithm::Impl {
    public:
        Impl(const MVNFileData& mvn);
        Impl(const vector<double>& lower,
             const vector<double>& upper,
             const vector< vector<double> >& sigma);
        bool checkAndPrepare();
        double integrand(const double x[]);
        double getUpperBound() const;
        uint32_t s;
        bool internalError;
        string errorString;
        Vectord lower;
        Vectord upper;
        Vectord xmu;
        Matrixd sigma;
        Matrixd L;
    };

    BotevAlgorithm::Impl::Impl(const MVNFileData& mvn)
    {
        internalError = false;
        s = mvn.s;
        lower.resize(s);
        upper.resize(s);
        sigma.resize(s, s);
        for (uint32_t i = 0; i < s; i++) {
            lower(i) = mvn.lower[i];
            upper(i) = mvn.upper[i];
        }
        for (uint32_t i = 0; i < s; i++) {
            for (uint32_t j = 0; j < s; j++) {
                sigma(i, j) = mvn.sigma[i * s + j];
            }
        }
    }

    BotevAlgorithm::Impl::Impl(const vector<double>& lower,
                               const vector<double>& upper,
                               const vector< vector<double> >& sigma)
    {
        internalError = false;
        s = lower.size();
        if (s != upper.size() || s != sigma.size() || s != sigma[0].size()) {
            internalError = true;
            errorString = "size mismatch";
            return;
        }
        this->lower.resize(s);
        this->upper.resize(s);
        this->sigma.resize(s, s);
        for (uint32_t i = 0; i < s; i++) {
            this->lower(i) = lower[i];
            this->upper(i) = upper[i];
        }
        for (uint32_t i = 0; i < s; i++) {
            for (uint32_t j = 0; j < s; j++) {
                this->sigma(i, j) = sigma[i][j];
            }
        }
    }

    BotevAlgorithm::~BotevAlgorithm()
    {
        delete impl;
    }

    BotevAlgorithm::BotevAlgorithm(uint32_t param_id, uint32_t s_dimR)
    {
        MVNFileData mvn(s_dimR);
        impl = new Impl(mvn);
        if (param_id >= mvn_parameter_name_data_size) {
            impl->internalError = true;
            impl->errorString = "param_id out of range";
            return;
        }
        int r = read_mvn(param_id, s_dimR, mvn);
        if (r != 0) {
            impl->internalError = true;
            impl->errorString = "can't read file";
        }
    }

    BotevAlgorithm::BotevAlgorithm(const vector<double>& lower,
                               const vector<double>& upper,
                               const vector< vector<double> >& sigma)
    {
        impl = new Impl(lower, upper, sigma);
    }

    double BotevAlgorithm::operator()(const double x[]) {
        return impl->integrand(x);
    }

    bool BotevAlgorithm::checkAndPrepare()
    {
        return impl->checkAndPrepare();
    }

    double BotevAlgorithm::Impl::getUpperBound() const
    {
        return upbnd(s, lower, upper, L, xmu);
    }

    double BotevAlgorithm::getUpperBound() const
    {
        return impl->getUpperBound();
    }

    bool BotevAlgorithm::Impl::checkAndPrepare()
    {
#if defined(DEBUG)
        cout << "in checkAndPrepare" << endl;
#endif
        if (internalError) {
            cout << errorString << endl;
            return false;
        }
#if defined(DEBUG)
        cout << "after internalError" << endl;
#endif
        L.resize(s, s);
        if (!choleskyDecomposition(s, sigma, L)) {
            return false;
        }
#if defined(DEBUG)
        cout << "after choleskeyDecomposition" << endl;
#endif
        int perm[s];
        cholPerm(s, sigma, lower, upper, perm, L);
#if defined(DEBUG)
        cout << "after cholPerm" << endl;
#endif
        rescale(s, lower, upper, L);
#if defined(DEBUG)
        cout << "after rescale" << endl;
#endif
        xmu.resize(2*s-2);
        xmu.clear();
        if (!nleq(s, lower, upper, L, xmu)) {
            return false;
        }
#if defined(DEBUG)
        cout << "after nleq" << endl;
#endif
        return true;
    }

/**
 * Integrand
 * use mvn.s
 *     mvn.lower
 *     mvn.upper
 *     xmu
 *     L
 * not use Jacobian
 *         grad
 */
    double BotevAlgorithm::Impl::integrand(const double x[])
    {
        Kahan ret;
        Kahan col;
        double Z[s];
        for (uint32_t i = 0; i < s; ++i) {
            Z[i] = 0.0;
        }
        for (uint32_t k = 0; k < s-1; ++k) {
            col.clear();
            for (uint32_t i = 0; i < k; ++i) {
                col.add(L(k, i) * Z[i]);
            }
            // mu(i) = xmu(s-1+i) for 0<=i<s-1
            double tl = lower(k) - xmu(s-1+k) - col.get();
            double tu = upper(k) - xmu(s-1+k) - col.get();
            Z[k] = xmu(s-1+k) + norminvp(x[k], tl, tu);
            ret.add(lnNpr(tl, tu) + xmu(s-1+k) * (0.5 * xmu(s-1+k) - Z[k]));
        }
        col.clear();
        for (uint32_t i = 0; i < s-1; ++i) {
            col.add(L(s-1, i) * Z[i]);
        }
        ret.add(lnNpr(lower(s-1) - col.get(), upper(s-1) - col.get()));
        return exp(ret.get());
    }


    uint32_t BotevAlgorithm::getParameterSize()
    {
        return mvn_parameter_name_data_size;
    }

    const string BotevAlgorithm::getParameterName(uint32_t index)
    {
        if (index < mvn_parameter_name_data_size) {
            return mvn_parameter_name_data[index].name;
        } else {
            return "";
        }
    }

    const string BotevAlgorithm::getParameterConstruction(uint32_t index)
    {
        if (index < mvn_parameter_name_data_size) {
            return mvn_parameter_name_data[index].construction;
        } else {
            return "";
        }
    }

}

#endif
