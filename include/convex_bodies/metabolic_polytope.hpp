// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef METABOLIC_POLYTOPE_HPP
#define METABOLIC_POLYTOPE_HPP

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <cmath>
#include <vector>

// This class describes a (generally not full dimensional) polytope defined by box bounds
// and equality constraints:  b_l <= x <= b_u, A_eq x = b_eq.
//
// This is a representation of a Metabolic Network. The variables correspond to reaction
// fluxes, b_l/b_u are the lower/upper flux bounds, and the equalities A_eq x = b_eq 
// encode the steady state condition S x = 0 (A_eq = S, b_eq = 0).
// @tparam Point Point type
// @tparam MT_Type the equality matrix type
template
<
    typename Point,
    typename MT_Type = Eigen::SparseMatrix<typename Point::FT, Eigen::RowMajor>
>
class MetabolicPolytope {
    public:
        typedef Point PointType;
        typedef MT_Type MT;
        typedef typename Point::FT NT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
        typedef Eigen::Triplet<NT> Triplet;
    private:
        unsigned d;  // the dimension
        MT A_eq;     // equality constraint matrix
        VT b_l, b_u; // lower & upper bounds on x
        VT b_eq;     // equality RHS, A_eq x = b_eq


    public:
        // Default constructor.
        MetabolicPolytope() = default;

        // @param d_ the dimension
        // @param A_eq_ the equality matrix, s.t. A_eq x = b_eq
        // @param b_l_ the lower bound vector
        // @param b_u_ the upper bound vector
        // @param b_eq_ the equality RHS vector
        MetabolicPolytope(unsigned d_, 
            MT const& A_eq_, 
            VT const& b_l_, 
            VT const& b_u_, 
            VT const& b_eq_
        ) : 
            d{d_}, A_eq{A_eq_}, b_l{b_l_}, b_u{b_u_}, b_eq{b_eq_}
        {}

        // Default copy constructor, copies all members.
        MetabolicPolytope(MetabolicPolytope const&) = default;

        // @return the dimension d
        unsigned getDimension() const { return this->d; }

        // @return the equality matrix A_eq
        MT const& getEqualities() const { return this->A_eq; }

        // @return the upper bound vector b_u
        VT const& getUpperBounds() const { return this->b_u; }

        // @return the lower bound vector b_l
        VT const& getLowerBounds() const { return this->b_l; }

        // @return the equality bound vector
        VT const& getEqualityBounds() const { return this->b_eq; }

        // @return true if the polytope has equality constraints
        bool hasEqualities() const { return A_eq.rows() > 0; }

        // @return the number of equality constraints (rows of A_eq)
        unsigned getNumEqualities() const { return (unsigned)A_eq.rows(); }

        // @return the number of finite bounds (finite entries in b_l and b_u)
        unsigned getNumFiniteBounds() const {
            unsigned count = 0;
            for (unsigned j = 0; j < d; ++j) {
                if (!std::isinf(b_l(j))) ++count;
                if (!std::isinf(b_u(j))) ++count;
            }
            return count;
        }
        
        // Builds a d dimensional cube with -1 <= xi <= 1,
        // for every dimension i.
        // @param d the dimension
        // @return the cube as a MetabolicPolytope
        static MetabolicPolytope cube(unsigned d) {
            VT b_l = -VT::Ones(d);
            VT b_u = VT::Ones(d);
            MT A_eq(0, d);
            VT b_eq = VT::Zero(0);
            return MetabolicPolytope(d, A_eq, b_l, b_u, b_eq);
        }

        // Builds a d dimensional simplex with x_i >= 0 for every
        // dimension i, and sum(x_i) = 1.
        // @param d the dimension
        // @return the simplex as a MetabolicPolytope
        static MetabolicPolytope simplex(unsigned d) {
            VT b_l = VT::Zero(d);
            VT b_u = VT::Ones(d);
            VT b_eq = VT::Ones(1);

            std::vector<Triplet> triplets;
            for (unsigned j = 0; j < d; ++j)
                triplets.push_back(Triplet(0, j, NT(1)));
                
            MT A_eq(1, d);
            A_eq.setFromTriplets(triplets.begin(), triplets.end());
            A_eq.makeCompressed();

            return MetabolicPolytope(d, A_eq, b_l, b_u, b_eq);
        }
};
#endif