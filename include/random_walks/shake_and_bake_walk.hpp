// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

/* EXPLANATION:

This is Running variant of Shake And Bake class of boundary sampling algorithms. 

[1] C. G. E. Boender, R. J. Caron, J. F. McDonald, A. H. G. Rinnooy Kan,  
    H. E. Romeijn, R. L. Smith, J. Telgen i A. C. F. Vorst,  
    *Shake-And-Bake Algorithms for Generating Uniform Points on the Boundary of Bounded Polyhedra*, 1991.  
    Available at: https://doi.org/10.1016/0166-218X(91)90006-7

*/

#ifndef RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <stdexcept> 

#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/correlation_matrices/corre_matrix.hpp"

struct ShakeAndBakeWalk
{

    template <typename Polytope, typename RandomNumberGenerator>
    struct Walk
    {
        using Point = typename Polytope::PointType;
        using VT = typename Polytope::VT;
        using NT = typename Point::FT;
        using MT = typename Polytope::MT;

        struct update_parameters 
        {
            int facet_prev = -1;
            NT  inner_vi_ak = NT(0);   
        };

        update_parameters _params;

        static constexpr NT kDefaultEpsilon = NT(1e-10);

        template <typename GenericPolytope>
        Walk(GenericPolytope& P,
            const Point& boundary_pt, 
            int facet_idx,    
            RandomNumberGenerator& rng,
            NT eps = kDefaultEpsilon)
            : epsilon_{eps}
        {
            P.normalize();
            initialize(P,boundary_pt, facet_idx, rng);
        }

        NT   get_epsilon() const noexcept { return epsilon_; }

        void apply(Polytope& P, unsigned int walk_len, RandomNumberGenerator& rng)
        {
            const NT eps = epsilon_; 

            for (unsigned step = 0; step < walk_len; ++step)
            {
                Point v = SBDirection<Point>::apply(P.dimension(),_A_row_k, rng);

                int facet_new;
                std::tie(_lambda_hit, facet_new) = P.line_positive_intersect(_p, v, _Ar, _Av, _lambda_hit, _params);

                if (!std::isfinite(_lambda_hit) || _lambda_hit <= eps  || facet_new < 0) 
                {
                    _lambda_hit = NT(0);
                    continue;
                }

                _p += _lambda_hit * v;
                _A_row_k = P.get_row(facet_new);
                _params.facet_prev  = facet_new;
            }
        }


        const Point& getCurrentPoint() const noexcept { return _p; }

    protected:

        void initialize(Polytope& P,
                        const Point& boundary_pt,
                        int facet_idx,
                        RandomNumberGenerator& rng)
        {
            int _dim = P.dimension();
            int _m = P.num_of_hyperplanes();
            VT b=P.get_vec();
            int active_facet;

            NT kFacetEps = epsilon_;

            // Checking if boundary point belongs to facet_idx
            _p = boundary_pt;
            VT ai = P.get_row(facet_idx);
            NT dist = std::abs(ai.dot(_p.getCoefficients()) - b.coeff(facet_idx));
            if (dist > kFacetEps)
            {
                active_facet = -1;
                for (int i = 0; i < _m; ++i) 
                {
                    VT ai = P.get_row(i);
                    NT dist = std::abs(ai.dot(_p.getCoefficients()) - b.coeff(i));
                    if (dist < kFacetEps) 
                    {
                        active_facet = i;
                        break;
                    }
                }
                if (active_facet < 0)
                {
                    throw std::runtime_error("Boundary point not on any facet!");
                }
            }
            active_facet = facet_idx;

            //Normal of active facet
            _A_row_k = P.get_row(active_facet);

            //Initializing Ar and Av 
            _Ar.setZero(_m);
            _Av.setZero(_m);
            _lambda_hit = NT(0);
            
            //Calculating first Ar and initializing lambda 
            _Ar.noalias() = P.get_mat() * _p.getCoefficients();
            _lambda_hit = NT(0);

            _A_row_k = P.get_row(active_facet);

            _params.facet_prev = active_facet;

        }

        NT epsilon_{kDefaultEpsilon};
        Point _p;
        VT _Ar;            
        VT _Av;            
        NT _lambda_hit;
        VT _A_row_k;
    };
};

#endif // RANDOM_WALKS_SHAKE_AND_BAKE_WALK_HPP