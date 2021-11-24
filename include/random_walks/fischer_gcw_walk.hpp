// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_FISCHER_GCW_WALK_HPP
#define RANDOM_WALKS_FISCHER_GCW_WALK_HPP


#include "sampling/sphere.hpp"
#include "sampling/sampling_segment.hpp"
#include <cmath>

// Random directions hit-and-run walk with uniform target distribution

struct FischerGCWalk
{
    struct parameters {};
    parameters param;

template
<
    typename Polytope,
    typename RandomNumberGenerator,
    typename LltDecomposition
>
struct Walk
{
    typedef typename Polytope::VT VT;
    typedef typename Polytope::MT MT;
    typedef typename Polytope::NT NT;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p, VT const& mu, NT const& k, unsigned int const& W, RandomNumberGenerator& rng)
    {
        _mu = mu;
        _sigma = MT::Identity(P.dimension(), P.dimension());
        _lltOfSigma = LltDecomposition(_sigma);
        _L_cov = _lltOfSigma.matrixL();
        _W = W;
        _k = k;
        initialize(P, p, k, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p, VT const& mu, MT const& sigma, NT const& k, unsigned int const& W,
         RandomNumberGenerator& rng, parameters const& params)
    {
        _mu = mu;
        _sigma = MT::Identity(P.dimension(), P.dimension());
        _lltOfSigma = LltDecomposition(_sigma);
        _L_cov = _lltOfSigma.matrixL();
        _W = W;
        _k = k;
        initialize(P, p, k, rng);
    }

    template
    <
        typename BallPolytope
    >
    inline void apply(BallPolytope const& P,
                      VT& p,   // a point to start
                      NT const &k,
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            GetGaussianDirectionTangentPlane<VT>::apply(p, _v, _L_cov, _sigma, rng);
            
            //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av, _lambda);

            _mu_p = _mu_p*cos(_lambda) + _mu_v*sin(_lambda);
            _mu_v = _mu.dot(_v);

            _lambda = sample_fischer_segment(bpair.second, bpair.first, _mu_p, _mu_v, _W, k, rng);

            p = (cos(_lambda) * p) + (sin(_lambda) * _v);
            //VT q = P.get_mat()*p - P.get_vec();
            //for (int i=0; i<P.num_of_hyperplanes(); i++)
            //{
            //    if (q(i)>NT(0))
            //   {
             //       std::cout<<"outside from sampling, q: "<<q(i)<<std::endl;
            //        exit(-1);
            //    }
            //}
        }
        //p = _p;
    }


    template
    <
        typename BallPolytope
    >
    inline void apply_ratio_esti(BallPolytope const& P,
                      VT& p,   // a point to start
                      NT const &k,
                      unsigned int const& walk_length,
                      RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            GetDirectionTangentPlane<VT>::apply(p, _v, rng);
            
            //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av, _lambda);

            _mu_p = _mu_p*cos(_lambda) + _mu_v*sin(_lambda);
            _mu_v = _mu.dot(_v);

            _lambda = sample_fischer_segment(bpair.second, bpair.first, _mu_p, _mu_v, _W, k, rng);

            p = (cos(_lambda) * p) + (sin(_lambda) * _v);
            //VT q = P.get_mat()*p - P.get_vec();
            //for (int i=0; i<P.num_of_hyperplanes(); i++)
            //{
            //    if (q(i)>NT(0))
            //   {
             //       std::cout<<"outside from sampling, q: "<<q(i)<<std::endl;
            //        exit(-1);
            //    }
            //}
        }
        //p = _p;
    }


    template
    <
        typename BallPolytope
    >
    inline void apply_with_check(BallPolytope const& P,
                                 VT& p,   // a point to start
                                 NT const &k,
                                 unsigned int const& walk_length,
                                 RandomNumberGenerator& rng)
    {
        for (auto j=0u; j<walk_length; ++j)
        {
            GetDirectionTangentPlane<VT>::apply(p, _v, rng);
            //_sigma_inv_v = _sigma_inv * _v;
            //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av, _lambda);

            _mu_p = _mu_p*cos(_lambda) + _mu_v*sin(_lambda);
            _mu_v = _mu.dot(_v);

            

            _res = sample_fischer_segment_with_check(bpair.second, bpair.first, _mu_p, _mu_v, _W, k, rng);
            _lambda = _res.first;
            _total_calls++;

            if (_res.second)
            {
                _is_out++;
                got_outside = true;
            }

            p = (cos(_lambda) * p) + (sin(_lambda) * _v);
        }
        //p = _p;
    }

    inline void set_sigma(MT const& sigma)
    {
        _sigma = sigma;
        _lltOfSigma = LltDecomposition(_sigma);
        _L_cov = _lltOfSigma.matrixL();
    }

    inline bool is_outside()
    {
        return got_outside;
    }

    inline bool ratio_outside()
    {
        return NT(_is_out) / NT(_total_calls);
    }

    template <typename BallPolytope>
    inline void initialize(BallPolytope const& P,
                           VT& p,
                           NT const &k,
                           RandomNumberGenerator &rng)
    {
        _lamdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _v.setZero(P.dimension());
        _is_out = 0;
        _total_calls = 0;
        got_outside = false;

        _mu_v = _mu.dot(_v);
        _mu_p = _mu.dot(p);

        //_sigma_p = _sigma * p;

        GetGaussianDirectionTangentPlane<VT>::apply(p, _v, _L_cov, _sigma, rng);

        //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
        std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av);
        _lambda = sample_fischer_segment(bpair.second, bpair.first, _mu_p, _mu_v, _W, k, rng);
        p = (cos(_lambda) * p) + (sin(_lambda) * _v);
    }

private :

    //Point _p;
    NT _lambda;
    VT _lamdas;
    VT _Av;
    VT _v;
    VT _mu;
    NT _mu_p;
    NT _mu_v;
    MT _sigma;
    NT _k;
    unsigned int _W, _is_out, _total_calls;
    bool got_outside = false;
    std::pair<NT, bool> _res;
    LltDecomposition _lltOfSigma;
    MT _L_cov; 
};

};


#endif // RANDOM_WALKS_UNIFORM_GCW_WALK_HPP
