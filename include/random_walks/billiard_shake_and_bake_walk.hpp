// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

/* EXPLANATION:

This is a variant of Running Shake and Bake algorithm implemented in 'shake_and_bake_walk.hpp' that belongs to the class of boundary sampling algorithms. 
It follows the steps as described in [1], but after step 1 (direction sampling) it does number of reflections defined by upper bound of reflection (nr). 
For each step, the number of reflections is sampled by using inverse exponential distribution, i.e. the number of reflections is (1-z)*nr. 

[1] C. G. E. Boender, R. J. Caron, J. F. McDonald, A. H. G. Rinnooy Kan,H. E. Romeijn, R. L. Smith, J. Telgen i A. C. F. Vorst,  
*Shake-And-Bake Algorithms for Generating Uniform Points on the Boundary of Bounded Polyhedra*, 1991.  
Available at: https://doi.org/10.1016/0166-218X(91)90006-7

*/

#ifndef RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP
#define RANDOM_WALKS_BILLIARD_SHAKE_AND_BAKE_WALK_HPP

#include <Eigen/Eigen>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <iostream>

#include "sampling/sphere.hpp"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/convex_body.h"
#include "random_walks/accelerated_billiard_walk_utils.hpp"
#include "random_walks/shake_and_bake_walk.hpp"
#include <boost/random/exponential_distribution.hpp>

struct BilliardShakeAndBakeWalk
{
    enum class ReflectionMode { Uniform, InverseExponential };

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(-1), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0), moved_dist(0.0) 
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
        double moved_dist; 
    };

    template
    <
            typename Polytope,
            typename RandomNumberGenerator
    >
    struct Walk: public ShakeAndBakeWalk::Walk<Polytope, RandomNumberGenerator>
    {
        using ShakeAndBake  = ShakeAndBakeWalk::Walk<Polytope, RandomNumberGenerator>;
        using Point = typename Polytope::PointType;
        using VT = typename Polytope::VT;
        using NT = typename Point::FT;
        using MT = typename Polytope::MT;
        using ShakeAndBake::initialize; // we initialize same as Shake and Bake 
        typedef typename Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DenseMT;

        static constexpr NT kDefaultEpsilon = NT(1e-10);
        static constexpr ReflectionMode kDefaultMode = ReflectionMode::InverseExponential;

        template <typename GenericPolytope>
        Walk(GenericPolytope &P, 
            Point const& p, 
            RandomNumberGenerator &rng,
            int facet_idx, 
            int nr, //upper bound for reflections
            NT eps = kDefaultEpsilon, 
            ReflectionMode mode = kDefaultMode): ShakeAndBakeWalk::template Walk<Polytope, RandomNumberGenerator>(P, p, facet_idx, rng, eps), mode_(mode)
        {
            if(!P.is_normalized()) 
            {
                P.normalize();
            }
            _params = update_parameters();

            // if not given, square root of dimension
            _nr = (nr > 0)? nr : static_cast<unsigned int>(std::ceil(std::sqrt(P.dimension())));

            _AA.noalias() = (DenseMT)(P.get_mat() * P.get_mat().transpose());

            //Initialize from Running Shake and Bake 
            initialize(P, p, facet_idx, rng);
        }

        NT get_epsilon() const noexcept { return this->epsilon_; }
        ReflectionMode get_mode() const noexcept { return this->mode_; }

        void apply(Polytope& P, unsigned int walk_len, RandomNumberGenerator& rng)
        {
            const NT eps = this->epsilon_;
            ReflectionMode mode = this->mode_;
            typename Point::Coeff b;

            for (unsigned int step = 0; step < walk_len; ++step)
            {
                _params.moved_dist = 0.0;
                unsigned int r;
                if (mode == ReflectionMode::Uniform)
                {
                    r = (_nr == 1) ? 1 : 1 + static_cast<unsigned int>(rng.sample_urdist() * _nr);
                }
                else // InverseExponential
                {
                    double z = rng.sample_trunc_expdist();
                    r = static_cast<unsigned int>(std::floor((1.0 - z) * _nr));
                }
 
                Point _v = SBDirection<Point>::apply(P.dimension(), this->_A_row_k, rng);
                auto pbair = P.line_first_positive_intersect(this->_p, _v, this->_Ar, this->_Av, _params);
                NT _lambda_prev = pbair.first;
                if (!std::isfinite(_lambda_prev) || _lambda_prev <= eps  || pbair.second < 0) 
                {
                    _lambda_prev = NT(0);
                    continue;
                }

                // from here same as accelerated billiard walk 

                this->_p += (_lambda_prev * _v);
   
                this->_A_row_k = P.get_row(pbair.second);
                _params.facet_prev = pbair.second;

                for (unsigned int k = 1; k < r; ++k) // from there we do reflections
                {
                    P.compute_reflection(_v, this->_p, _params);
                    pbair = P.line_positive_intersect(this->_p, _v, this->_Ar, this->_Av, _lambda_prev, _AA, _params);


                    _lambda_prev = pbair.first;
                    if (!std::isfinite(_lambda_prev) || _lambda_prev <= eps  || pbair.second < 0) 
                    {
                        _lambda_prev = NT(0);
                        continue;
                    }
                
                    _params.moved_dist += _lambda_prev;

                    this->_p += _lambda_prev * _v;

                    this->_A_row_k = P.get_row(pbair.second);       
                    _params.facet_prev = pbair.second;
                }
            }
        }

        const Point& getCurrentPoint() const noexcept { return this->_p; }
        
    private:
    
        DenseMT _AA;
        update_parameters _params;
        BoundaryOracleHeap<NT> _distances_set;
        int _nr; 
        ReflectionMode mode_;
    };

};

#endif