// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GAUSSIAN_GCW_WALK_HPP
#define RANDOM_WALKS_GAUSSIAN_GCW_WALK_HPP


#include "sampling/sphere.hpp"
#include "sampling/sampling_segment.hpp"
#include <cmath>

// Random directions hit-and-run walk with uniform target distribution

struct GaussianGCWalk
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
    Walk(GenericPolytope const& P, VT& p, VT const& mu, MT const& sigma, NT const& k, unsigned int const& W, RandomNumberGenerator& rng)
    {
        _mu = mu;
        _sigma_inv = sigma.inverse();
        _sigma = sigma;
        _lltOfSigma = LltDecomposition(sigma);
        _L_cov = _lltOfSigma.matrixL();
        _sigma_inv_mu = _sigma_inv * mu;
        _W = W;
        _k = k;
        initialize(P, p, k, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, VT& p, VT const& mu, MT const& sigma, NT const& k, unsigned int const& W,
         RandomNumberGenerator& rng, parameters const& params)
    {
        _mu = mu;
        _sigma_inv = sigma.inverse();
        _sigma = sigma;
        _lltOfSigma = LltDecomposition(sigma);
        _L_cov = _lltOfSigma.matrixL();
        _sigma_inv_mu = _sigma_inv * mu;
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
            _sigma_inv_v = _sigma_inv * _v;
            std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av, _lambda);

            a = a*cos(_lambda) + c*cos(_lambda)*sin(_lambda) + b*sin(_lambda)*sin(_lambda);
            b = _v.dot(_sigma_inv_v);
            c = NT(2) * (p.dot(_sigma_inv_v));
            d = -NT(2) * (_v.dot(_sigma_inv_mu));
            e = -NT(2) * (p.dot(_sigma_inv_mu)); // todo: optize it

            _lambda = sample_sigma_gaussian_segment(bpair.second, bpair.first, a, b, c, d, e, _W, k, rng);

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
            GetGaussianDirectionTangentPlane<VT>::apply(p, _v, _L_cov, _sigma, rng);
            _sigma_inv_v = _sigma_inv * _v;
            std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
            std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av, _lambda);

            a = a*cos(_lambda) + c*cos(_lambda)*sin(_lambda) + b*sin(_lambda)*sin(_lambda);
            b = _v.dot(_sigma_inv_v);
            c = NT(2) * (p.dot(_sigma_inv_v));
            d = -NT(2) * (_v.dot(_sigma_inv_mu));
            e = -NT(2) * (p.dot(_sigma_inv_mu)); // todo: optize it

            _lambda = sample_sigma_gaussian_segment(bpair.second, bpair.first, a, b, c, d, e, _W, k, rng);

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

    inline bool is_inside()
    {
        return got_outside;
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

        _sigma_inv_v = _sigma_inv * _v;
        //_sigma_inv_p = _sigma_inv * p;
        a = p.dot(_sigma_inv * p);
        b = _v.dot(_sigma_inv_v);
        c = NT(2) * (p.dot(_sigma_inv_v));
        d = -NT(2) * (_v.dot(_sigma_inv_mu));
        e = -NT(2) * (p.dot(_sigma_inv_mu));

        //_sigma_p = _sigma * p;

        GetGaussianDirectionTangentPlane<VT>::apply(p, _v, _L_cov, _sigma, rng);

        //std::cout<<"p'v = "<<p.dot(_v)<<", v.norm() = "<<_v.norm()<<std::endl;
        std::pair<NT, NT> bpair = P.gc_intersect(p, _v, _lamdas, _Av);
        _lambda = sample_sigma_gaussian_segment(bpair.second, bpair.first, a, b, c, d, e, _W, k, rng);
        p = (cos(_lambda) * p) + (sin(_lambda) * _v);
    }

private :

    //Point _p;
    NT _lambda;
    VT _lamdas;
    VT _Av;
    VT _v;
    VT _mu;
    VT _sigma_inv_mu;
    VT _sigma_inv_v;
    //VT _sigma_inv_p;
    MT _sigma_inv;
    MT _sigma;
    NT a, b, c, d, e, _k;
    unsigned int _W;
    bool got_outside = false;
    std::pair<NT, bool> res;
    LltDecomposition _lltOfSigma;
    MT _L_cov; 
};

};


#endif // RANDOM_WALKS_UNIFORM_GCW_WALK_HPP
