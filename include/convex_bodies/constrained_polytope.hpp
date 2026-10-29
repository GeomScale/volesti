// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CONSTRAINED_POLYTOPE_HPP
#define CONSTRAINED_POLYTOPE_HPP

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>

// This class describes a (generally not full dimensional) polytope defined by equality 
// constraints, general inequality constraints and bounds on variables:
//
//      A_eq x = b_eq
//      A_in x <= b_in
//      lb <= x <= ub
//
// @tparam Point Point type
template <typename Point>
class ConstrainedPolytope {
    public:
        typedef Point PointType;
        typedef typename Point::FT NT;
        typedef Eigen::SparseMatrix<NT, Eigen::RowMajor> MT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
        typedef Eigen::Triplet<NT> Triplet;

    private:
        unsigned d;  // the ambient dimension
        MT A_eq;     // equality constraint matrix
        VT b_eq;     // equality rhs
        MT A_in;     // inequality matrix
        VT b_in;     // inequality rhs
        VT lb;       // lower bounds on x
        VT ub;       // upper bound on x

    public:
        // Default constructor.
        ConstrainedPolytope() = default;

        // Builds a polytope with equalities, general inequalities and box bounds.
        // @param d_ the ambient dimension
        // @param A_eq_ the equality matrix, s.t. A_eq x = b_eq
        // @param b_eq_ the equality rhs
        // @param A_in_ the inequality matrix
        // @param b_in_ the inequality rhs
        // @param lb_ the lower bounds on x
        // @param ub_ the upper bounds on x
        ConstrainedPolytope(unsigned d_, 
                          MT const& A_eq_, VT const& b_eq_, 
                          MT const& A_in_,  VT const& b_in_,
                          VT const& lb_, VT const& ub_
        ) : 
            d{d_}, A_eq{A_eq_}, b_eq{b_eq_}, A_in{A_in_}, b_in{b_in_}, lb{lb_}, ub{ub_}
        {}

        // Builds a polytope with equalities, general inequalities and no box bounds.
        // @param d_ the ambient dimension
        // @param A_eq_ the equality matrix, s.t. A_eq x = b_eq
        // @param b_eq_ the equality rhs
        // @param A_in_ the inequality matrix
        // @param b_in_ the inequality rhs
        ConstrainedPolytope(unsigned d_, 
                          MT const& A_eq_, VT const& b_eq_, 
                          MT const& A_in_,  VT const& b_in_
        ) : 
            ConstrainedPolytope(d_, A_eq_, b_eq_, A_in_, b_in_,
                                VT::Constant(d_, -std::numeric_limits<NT>::infinity()),
                                VT::Constant(d_, std::numeric_limits<NT>::infinity()))
        {}

        // Default copy constructor, copies all members.
        ConstrainedPolytope(ConstrainedPolytope const&) = default;

        // @return the dimension d
        unsigned getDimension() const { return this->d; }

        // @return the equality matrix A_eq
        MT const& getEqualities() const { return this->A_eq; }

        // @return the equality rhs
        VT const& getEqualityRHS() const { return this->b_eq; }

        // @return the inequality matrix A_in
        MT const& getInequalities() const { return this->A_in; }

        // @return the inequality rhs
        VT const& getInequalityRHS() const { return this->b_in; }

        // @return the lower bounds on x
        VT const& getLowerBounds() const { return this->lb; }

        // @return the upper bounds on x
        VT const& getUpperBounds() const { return this->ub; }

        // @return true if the polytope has equality constraints
        bool hasEqualities() const { return A_eq.rows() > 0; }

        // @return true if the polytope has inequality constraints
        bool hasInequalities() const { return A_in.rows() > 0; }

        // @return the number of equality constraints (rows of A_eq)
        unsigned getNumEqualities() const { return (unsigned)A_eq.rows(); }

        // @return the number of inequality constraints (rows of A_in)
        unsigned getNumInequalities() const { return (unsigned)A_in.rows(); }

        // @return the number of finite entries in lb
        unsigned getNumFiniteLowerBounds() const {
            unsigned n = 0;
            for (Eigen::Index j = 0; j < lb.size(); ++j)
                if (std::isfinite((double)lb(j))) ++n;

            return n;
        }

        // @return the number of finite entries in ub
        unsigned getNumFiniteUpperBounds() const {
            unsigned n = 0;
            for (Eigen::Index j = 0; j < ub.size(); ++j)
                if (std::isfinite((double)ub(j))) ++n;

            return n;
        }

        // @return the number of finite entries in lb, ub
        unsigned getNumFiniteBounds() const { return this->getNumFiniteLowerBounds()+this->getNumFiniteUpperBounds(); }
};
#endif