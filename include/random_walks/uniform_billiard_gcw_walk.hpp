// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BILLIARD_GCW_WALK_HPP
#define RANDOM_WALKS_UNIFORM_BILLIARD_GCW_WALK_HPP


#include "sampling/sphere.hpp"
#include <cmath>

// Billiard walk for uniform distribution

struct BilliardGCWalk
{
    BilliardGCWalk(double L)
            :   param(L, true)
    {}

    BilliardGCWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

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
    Walk(GenericPolytope const& P, VT &p, NT const& L, RandomNumberGenerator &rng)
    {
        _Len = L;
        _Identity_mat = MT::Identity(P.dimension(), P.dimension());
        _p_p_tr.setZero(P.dimension(), P.dimension());
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT &p, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.m_L;
        _Identity_mat = MT::Identity(P.dimension(), P.dimension());
        _p_p_tr.setZero(P.dimension(), P.dimension());
        initialize(P, p, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      VT& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T = rng.sample_urdist() * _Len;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;
            _p_p_tr = _Identity_mat - p * p.transpose();
            GetDirectionTangentPlaneOpt<VT, MT>::apply(p, _v, _p_p_tr, rng);
            _p0 = p;
            int it = 0;
            while (it < 50000*n)
            {
                _pbpair = P.gc_intersect_positive(p, _v, _lambdas, _Av, _lambda);
                if (T <= _pbpair.first || _pbpair.second < 0) {
                    p = (cos(T) * p) + (sin(T) * _v);
                    _lambda = T;
                    break;
                }
                _lambda = dl * _pbpair.first;
                _vi = (-sin(_lambda) * p) + (cos(_lambda) * _v);
                p = (cos(_lambda) * p) + (sin(_lambda) * _v);
                T -= _lambda;

                //if (P.is_in(p)==0){
                //    std::cout<<"point out"<<std::endl;
                //    exit(-1);
                //}
                
                _p_p_tr = _Identity_mat - p * p.transpose();
                P.compute_reflection(_vi, p, _p_p_tr, _pbpair.second);
                _v = _vi;

                _v = (_p_p_tr * _v).eval();
                _v *= (NT(1)/_v.norm());
                it++;
            }
            //if (P.is_in(p)==0){
            //    std::cout<<"point out"<<std::endl;
            //    exit(-1);
            //}
            if (it == 50000*n){
                p = _p0;
            }
        }
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }



    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope const& P,
                           VT &p,
                           RandomNumberGenerator &rng)
    {
        //std::cout<<"initializing"<<std::endl;
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _v.setZero(P.dimension());

        _p_p_tr = _Identity_mat - p * p.transpose();
        GetDirectionTangentPlaneOpt<VT, MT>::apply(p, _v, _p_p_tr, rng);

        NT T = rng.sample_urdist() * _Len;
        int it = 0;

        _pbpair = P.gc_intersect_positive(p, _v, _lambdas, _Av);
        if (T <= _pbpair.first || _pbpair.second < 0) {
            p = (cos(T) * p) + (sin(T) * _v);
            _lambda = T;
            return;
        }
        _lambda = dl * _pbpair.first;
        _vi = (-sin(_lambda) * p) + (cos(_lambda) * _v);
        p = (cos(_lambda) * p) + (sin(_lambda) * _v);
        T -= _lambda;

        //std::cout<<"initialized [1]"<<std::endl;
        _p_p_tr = _Identity_mat - p * p.transpose();
        //std::cout<<"initialized [2]"<<std::endl;
        P.compute_reflection(_vi, p, _p_p_tr, _pbpair.second);
        _v = _vi;
        //std::cout<<"initialized [3]"<<std::endl;
        //if (P.is_in(p)==0){
        //    std::cout<<"point out"<<std::endl;
         //   exit(-1);
        //}

        while (it <= 50000*n)
        {
            _v = (_p_p_tr * _v).eval();
            _v *= (NT(1)/_v.norm());

            _pbpair = P.gc_intersect_positive(p, _v, _lambdas, _Av, _lambda);

            if (T <= _pbpair.first || _pbpair.second < 0) {
                p = (cos(T) * p) + (sin(T) * _v);
                _lambda = T;
                break;
            }else if (it == 50000*n) {
                _lambda = rng.sample_urdist() * _pbpair.first;
                p = (cos(_lambda) * p) + (sin(_lambda) * _v);
                break;
            }
            _lambda = dl * _pbpair.first;
            _vi = (-sin(_lambda) * p) + (cos(_lambda) * _v);
            p = (cos(_lambda) * p) + (sin(_lambda) * _v);
            T -= _lambda;

            //if (P.is_in(p)==0){
            //    std::cout<<"point out"<<std::endl;
            //    exit(-1);
            //}

            _p_p_tr = _Identity_mat - p * p.transpose();
            P.compute_reflection(_vi, p, _p_p_tr, _pbpair.second);
            _v = _vi;
            it++;
        }
        //if (P.is_in(p)==0){
        //    std::cout<<"point out"<<std::endl;
        //    exit(-1);
        //}
    }

private :

    NT _Len;
    VT _p0;
    VT _v;
    VT _vi;
    NT _lambda;
    MT _Identity_mat;
    MT _p_p_tr;
    std::pair<NT, int> _pbpair;
    VT _lambdas;
    VT _Av;
};

};








#endif // RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP
