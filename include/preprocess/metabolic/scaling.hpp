// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_SCALING_HPP
#define METABOLIC_SCALING_HPP

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <limits>
#include "convex_bodies/metabolic_polytope.hpp"

// A scaling of the polytope. Both vectors are strictly positive and
// are divisors of the coefficients of A_eq, so that
// 
// A_eq*(i,j) = A_eq(i,j)/(row(i)*col(j)),
// b_eq*(i)   = b_eq(i)/row(i),
// x*(j)      = col(j)*x(j),
// b_l*(j)    = col(j)*b_l(j),
// b_u*(j)    = col(j)*b_u(j)
// @tparam Point the point type of the polytope
template <typename Point>
struct Scaling {
    typedef typename MetabolicPolytope<Point>::VT VT;
    VT col; // Scaling factors of the reactions.
    VT row; // Scaling factors of the metabolites.
};

// Rounds a scaling factor to the nearest power of two, in an effort
// to minimize floating point inaccuracies when multiplying.
// @param v the factor
// @return the nearest power of two, or 1 if v is not positive or finite
inline double round_pow2(double v) {
    double u = (!(v > 0.0) || std::isinf(v)) ? 1.0 : std::exp2(std::round(std::log2(v)));
    return u;
}

// Scales every reaction by its largest finite bound in absolute value,
// which sends that bound to +-1, and every metabolite row by its largest
// coefficieint once the reaction factors are in place, which sends that coeffecient
// to +-1.
struct MaxBoundScaling {
    // Sets the scaling factors for the polytope.
    // @tparam Point the point type of the polytope
    // @param P the polytope
    // @param s the scaling factors, sized to d and m
    template<typename Point>
    void operator()(MetabolicPolytope<Point> const& P, Scaling<Point> & s) const {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::VT VT;
        typedef typename MetabolicPolytope<Point>::NT NT;

        MT const& A_eq = P.getEqualities();
        VT const& b_u = P.getUpperBounds();
        VT const& b_l = P.getLowerBounds();
        unsigned d = P.getDimension();
        unsigned m = (unsigned)A_eq.rows();

        for (unsigned j = 0; j < d; ++j) {
            double a = 0.0;
            double b_uj = (double)b_u(j);
            double b_lj = (double)b_l(j);
            a = std::isfinite(b_uj) ? std::max(a, std::abs(b_uj)) : a;
            a = std::isfinite(b_lj) ? std::max(a, std::abs(b_lj)) : a;
            s.col(j) = (NT)round_pow2(a > 0.0 ? 1.0/a : 1.0);
        }

        for (unsigned i = 0; i < m; ++i) {
            double a = 0.0;
            for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                a = std::max(a, std::abs((double)it.value()/(double)s.col((unsigned)it.col())));
            }
            s.row(i) = (NT)round_pow2(a > 0.0 ? a : 1.0);
        }
    }
};

// Geometric mean scaling of A_eq, following the gmscale
struct GMScaling {
    // Number of alternating column/row passes.
    unsigned passes = 5;

    // Convergence tolerance. If a pass fails to improve the max/min ratio by this factor
    // stops the iteration.
    double scltol = 0.9;

    // Damping of the smallest magnitude in a row or column.
    double damp = 1e-4;

    // Scale factors below this are treated as degenerate and reset to one.
    double tol = 1e-12;

    // Sets the scaling factors for the polytope.
    // @tparam Point the point type of the polytope
    // @param P the polytope
    // @param s the scaling factors, sized to d and m
    template<typename Point>
    void operator()(MetabolicPolytope<Point> const& P, Scaling<Point> & s) const {
        typedef typename MetabolicPolytope<Point>::MT MT;
        typedef typename MetabolicPolytope<Point>::NT NT;

        MT const& A_eq = P.getEqualities();
        unsigned d = P.getDimension();
        unsigned m = (unsigned)A_eq.rows();

        std::vector<double> col(d, 1.0), row(m, 1.0), max_c, min_c;

        const double INF = std::numeric_limits<double>::infinity();
        const double EPS = 2.2204e-16;
        double aratio = 1e50;

        for (unsigned k = 0; k < passes; ++k) {
            max_c.assign(d, 0.0);
            min_c.assign(d, INF);

            for (unsigned i = 0; i < m; ++i) {
                for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                    double v = std::abs((double)it.value())/(double)row[i];
                    if (!v) continue;
                    unsigned j = (unsigned)it.col();
                    max_c[j] = std::max(max_c[j], v);
                    min_c[j] = std::min(min_c[j], v); 
                }
            }

            double sratio = 0.0;
            for (unsigned j = 0; j < d; ++j)
                if (max_c[j] > 0.0) sratio = std::max(sratio, max_c[j]*(1.0/(1.0/min_c[j]+EPS)));

            if (k > 0) {
                for (unsigned j = 0; j < d; ++j) {
                    if (max_c[j] <= 0.0) continue;
                    col[j] = std::sqrt(std::max(min_c[j], damp*max_c[j])*max_c[j]);
                }
            }

            if (k >= 2 && sratio >= aratio*scltol) 
                break;

            aratio = sratio;

            for (unsigned j = 0; j < d; ++j) 
                if (col[j] < tol) col[j] = 1.0;

            max_c.assign(m, 0.0);
            min_c.assign(m, INF);

            for (unsigned i = 0; i < m; ++i) {
                for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                    double v = std::abs((double)it.value())/col[(unsigned)it.col()];
                    if (v <= 0.0) continue;
                    min_c[i] = std::min(min_c[i], v);
                    max_c[i] = std::max(max_c[i], v);
                }
            }

            for (unsigned i = 0; i < m; ++i) {
                if (max_c[i] <= 0.0) continue;
                row[i] = std::sqrt(std::max(min_c[i], damp*max_c[i])*max_c[i]);
            }
        }

        for (unsigned i = 0; i < m; ++i)
            if (row[i] == 0.0) row[i] = 1.0;

        max_c.assign(d, 0.0);
        for (unsigned i = 0; i < m; ++i) {
            for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
                unsigned j = (unsigned)it.col();
                max_c[j] = std::max(max_c[j], std::abs((double)it.value())/row[i]);
            }
        }

        for (unsigned j = 0; j < d; ++j)
            col[j] = max_c[j] > 0.0 ? max_c[j] : 1.0;

        for (unsigned j = 0; j < d; ++j) s.col(j) = (NT)round_pow2(col[j]);
        for (unsigned i = 0; i < m; ++i) s.row(i) = (NT)round_pow2(row[i]);
    }
};


// Leaves the metabolic polytope untouched, useful for debugging purposes.
struct NoScaling {
    // Sets the scaling factors for the polytope.
    // @tparam Point the point type of the polytope
    // @param P the polytope
    // @param s the scaling factors, sized to d and m
    template <typename Point>
    void operator()(MetabolicPolytope<Point> const& P, Scaling<Point> & s) const {
        s.col.setOnes();
        s.row.setOnes();
    }
};

// Reconstructs the polytope under a scaling. The dimension, the number of
// equalities, the sparsity pattern and every index are left untouched, only
// numerical values move.
// @tparam Point the point type of the polytope
// @param P the metabolic polytope
// @param s the scaling
// @param invert if true the scaling is undone instead of being applied
// @return the scaled polytope
template <typename Point>
MetabolicPolytope<Point> rescale(MetabolicPolytope<Point> const& P,
                                    Scaling<Point> const& s,
                                    bool invert)
{
    typedef typename MetabolicPolytope<Point>::MT MT;
    typedef typename MetabolicPolytope<Point>::VT VT;
    typedef typename MetabolicPolytope<Point>::NT NT;

    MT A_eq = P.getEqualities();
    VT b_eq = P.getEqualityBounds();
    VT b_u = P.getUpperBounds();
    VT b_l = P.getLowerBounds();
    unsigned d = P.getDimension();
    unsigned m = (unsigned)A_eq.rows();

    for (unsigned i = 0; i < m; ++i) {
        for (typename MT::InnerIterator it(A_eq, i); it; ++it) {
            NT k = s.row(i)*s.col((unsigned)it.col());
            if (invert) {
                it.valueRef() *= k;
            } else {
                it.valueRef() /= k;
            }
        }
    }

    for (unsigned i = 0; i < m; ++i)
        b_eq(i) = invert ? b_eq(i)*s.row(i) : b_eq(i)/s.row(i);

    for (unsigned j = 0; j < d; ++j) {
        b_l(j) = invert ? b_l(j)/s.col(j) : b_l(j)*s.col(j);
        b_u(j) = invert ? b_u(j)/s.col(j) : b_u(j)*s.col(j);
    }

    return MetabolicPolytope<Point>(d, A_eq, b_l, b_u, b_eq);
}

// Maps the polytope into the scaled coordinates, using the given scaling
// policy to choose the factors.
//
// A policy is any callable with the signature
//
// void(MetabolicPolytope<Point> const& P, Scaling<Point> &s)
//
// The policy is expected to write strictly positive finite values into the entries
// of s
// @tparam Point the point type of the factor of the polytope
// @param s set to the scaling that will be applied
// @param policy the function that fills s
// @return the scaled polytope
template <typename Point, typename ScalingPolicy = MaxBoundScaling>
MetabolicPolytope<Point> scale(MetabolicPolytope<Point> const& P,
                                Scaling<Point> & s,
                                ScalingPolicy policy = ScalingPolicy{})
{
    static_assert(
        std::is_invocable_v<ScalingPolicy&, MetabolicPolytope<Point> const&, Scaling<Point>&>,
        "A ScalingPolicy must be a callable as void(MetabolicPolytope<Point> const&, Scaling<Point>&)");

    typedef typename MetabolicPolytope<Point>::VT VT;

    unsigned d = P.getDimension();
    unsigned m = (unsigned)P.getEqualities().rows();

    s.col = VT(d);
    s.row = VT(m);

    policy(P, s);

    return rescale(P, s, false);
}

// Maps a point of the original polytope into the scaled coordinates or
// maps a point of the scaled polytope back into the original coordinates
// depending on the flag is_original.
// @tparam Point the point type of the polytope
// @tparam XT the point type of x
// @param p the point
// @param s the scaling
// @param is_original true if the point is from the original polytope
// @return the scaled point
template <typename Point, typename XT>
XT scale_point(XT const& x, Scaling<Point> const& s, bool is_original) {
    XT y = x;
    for (unsigned j = 0; j < (unsigned)y.size(); ++j)
        y(j) = is_original ? x(j)*(typename XT::Scalar)s.col(j) 
                            : x(j)/(typename XT::Scalar)s.col(j);
    
    return y;
}

// This is a helper function that extends the metabolite factors with ones so that
// scaling stays valid after simplification appends equality rows for the fixed dimensions.
// @tparam Point the point type of the polytope
// @param the scaling to extend
// @param row_count the new number of equalities
template <typename Point>
void pad_rows(Scaling<Point> & s, unsigned row_count) {
    typedef typename MetabolicPolytope<Point>::VT VT;

    unsigned m = (unsigned)s.row.size();

    if (row_count <= m) // Guards the case in which simplification didn't add rows
        return;

    VT temp(row_count);
    temp.head(m) = s.row;
    temp.tail(row_count-m).setOnes();
    s.row = temp;
};  
#endif