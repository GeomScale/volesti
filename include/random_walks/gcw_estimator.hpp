// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_GCW_ESTIMATOR_HPP
#define RANDOM_WALKS_UNIFORM_GCW_ESTIMATOR_HPP


#include "sampling/sphere.hpp"
#include <cmath>

// Random directions hit-and-run walk with uniform target distribution

struct GCWEstimator
{
    struct parameters {};
    parameters param;

template
<
    typename Polytope,
    typename PointList,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::NT NT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p, RandomNumberGenerator& rng)
    {
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p,
         RandomNumberGenerator& rng, parameters const& params)
    {
        initialize(P, p, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      VT& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            VT v = GetDirectionTangentPlane<VT>::apply(p, rng);
            std::pair<NT, NT> bpair = P.gc_intersect(p, v, _lamdas, _Av,
                                                       _lambda);
            _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                    + bpair.second;
            p = (cos(_lambda) * p) + (sin(_lambda) * v);
        }
        //p = _p;
    }


    template
    <
        typename BallPolytope
    >
    inline PointList estimate(BallPolytope const& P1,
                              BallPolytope const& P2,
                              VT& p,   // a point to start
                              unsigned int const& walk_length,
                              NT const& error,
                              NT &val,
                              const unsigned int &W,
                              const unsigned int &Ntot,
                              const NT &ratio,
                              RandomNumberGenerator& rng)
    {
        int n = P1.dimension(), min_index = W-1, max_index = W-1, index = 0, iter = 1;
        PointList list_of_samples;
        std::vector<NT> last_W(W,0);
        size_t totCount = Ntot, countIn = Ntot * ratio;
        typename std::vector<NT>::iterator minmaxIt;
        NT min_val = std::numeric_limits<NT>::lowest(), max_val = std::numeric_limits<NT>::max();
        VT v(P1.dimension());
        bool verbose = true;

        P2.is_in_optimized(p, _lamdas, _Av, _lambda); //preprocessing

        while(iter <= MAX_ITER_ESTI){
            iter++;

            v = GetDirectionTangentPlane<VT>::apply(p, rng);
            std::pair<NT, NT> bpair = P1.gc_intersect_optimized(p, v, _lamdas, _Av, _lambda);
            _lambda = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
            p = (cos(_lambda) * p) + (sin(_lambda) * v);

            if(P2.is_in_optimized(p, _lamdas, _Av, _lambda)==-1){
                countIn++;

                if (store_points)
                {
                    list_of_samples.push_back(p);
                }
            }

            totCount++;
            val = NT(countIn) / NT(totCount);
            last_W[index] = val;

            if(val<=min_val){
                min_val = val;
                min_index = index;
            } else if(min_index==index){
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if(val>=max_val){
                max_val = val;
                max_index = index;
            } else if(max_index==index){
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            if( (max_val-min_val)/max_val<=error/2.0 ){
                if (verbose) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
                return list_of_samples;
            }

            index = index%W+1;
            if(index==W) index=0;

        }
        return list_of_samples;
    }


    inline void activate_storing()
    {
        store_points = true;
    }

    inline void deactivate_storing()
    {
        store_points = false;
    }

private :

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           VT& p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());

        VT v = GetDirectionTangentPlane<VT>::apply(p, rng);
        std::pair<NT, NT> bpair = P.gc_intersect(p, v, _lamdas, _Av);
        _lambda = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
        p = (cos(_lambda) * p) + (sin(_lambda) * v);
    }

    //Point _p;
    NT _lambda;
    VT _lamdas;
    VT _Av;
    bool store_points = false;
    const unsigned int MAX_ITER_ESTI = 80000000;
};

};


#endif // RANDOM_WALKS_UNIFORM_GCW_WALK_HPP
