// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_GCW_WALK_HPP
#define RANDOM_WALKS_UNIFORM_GCW_WALK_HPP


#include "sampling/sphere.hpp"
#include <cmath>

// Random directions hit-and-run walk with uniform target distribution

struct GCWalk
{
    struct parameters {};
    parameters param;

template
<
    typename Polytope,
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
            GetDirectionTangentPlane<VT>::apply(p, _v, rng);
            //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av,
                                                       _lambda);
            _lambda = rng.sample_urdist() * (bpair.first - bpair.second)
                    + bpair.second;
            p = (cos(_lambda) * p) + (sin(_lambda) * _v);
            //VT q = P.get_mat()*p - P.get_vec();
            //for (int i=0; i<P.num_of_hyperplanes(); i++)
            //{
            //    if (q(i)>NT(0))
            //   {
            //        std::cout<<"outside from sampling, q: "<<q(i)<<std::endl;
            //        exit(-1);
            //    }
            //}
        }
        //p = _p;
    }



    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           VT& p,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _v.setZero(P.dimension());

        GetDirectionTangentPlane<VT>::apply(p, _v, rng);
        //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
        std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av);
        _lambda = rng.sample_urdist() * (bpair.first - bpair.second) + bpair.second;
        p = (cos(_lambda) * p) + (sin(_lambda) * _v);
    }

private :

    //Point _p;
    NT _lambda;
    VT _lamdas;
    VT _Av;
    VT _v;
};

};


#endif // RANDOM_WALKS_UNIFORM_GCW_WALK_HPP
