// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP
#define RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP

#include <Eigen/Eigen>

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/random_point_generators.hpp"
#include "sampling/sphere.hpp"
#include "volume/sampling_policies.hpp"

template <typename GenericPolytope>
struct compute_diameter
{
    template <typename NT>
    static NT compute(GenericPolytope) {}
};


template <typename Point>
struct compute_diameter<HPolytope<Point>>
{
template <typename NT>
static NT compute(HPolytope<Point> const& P)
{
    NT diameter = NT(4) * std::sqrt(NT(P.dimension())) * P.InnerBall().second;
    return diameter;
}
};

template <typename Point>
struct compute_diameter<VPolytope<Point>>
{
template <typename NT>
static NT compute(VPolytope<Point> const& P)
{
    typedef typename VPolytope<Point>::MT MT;
    NT diameter = NT(0), diam_iter;
    MT V = P.get_mat();
    for (int i = 0; i < V.rows(); ++i) {
        for (int j = 0; j < V.rows(); ++j) {
            if (i != j) {
                diam_iter = (V.row(i) - V.row(j)).norm();
                if (diam_iter > diameter) diameter = diam_iter;
            }
        }
    }
    return diameter;
}
};

template <typename Point>
struct compute_diameter<Zonotope<Point>>
{
template <typename NT>
static NT compute(Zonotope<Point> const& P)
{
    typedef typename Zonotope<Point>::MT MT;
    typedef typename Zonotope<Point>::VT VT;

    MT V = P.get_mat();
    int k = V.rows(), max_index = -1;
    MT D = V.transpose() * V;
    D = (D + D.transpose()) / 2.0;
    Eigen::SelfAdjointEigenSolver <MT> es(D);
    MT D2 = es.eigenvalues().asDiagonal(), Q = es.eigenvectors();

    NT max_eig = NT(0);
    for (int i = 0; i < P.dimension(); ++i) {
        if (es.eigenvalues()[i] > max_eig) {
            max_eig = es.eigenvalues()[i];
            max_index = i;
        }
    }

    VT max_eigvec = -1.0 * Q.col(max_index);
    VT obj_fun = max_eigvec.transpose() * V.transpose(), x0(k);

    for (int j = 0; j < k; ++j) x0(j) = (obj_fun(j) < 0.0) ? -1.0 : 1.0;

    NT diameter = NT(2) * (V.transpose() * x0).norm();
    return diameter;
}
};

template <typename Point, typename RandomNumberGenerator>
struct compute_diameter<IntersectionOfVpoly<VPolytope<Point>, RandomNumberGenerator>>
{
template <typename NT>
static NT compute(IntersectionOfVpoly<VPolytope<Point>, RandomNumberGenerator> const& P)
{
    NT diameter = NT(2) * NT(P.dimension()) * P.InnerBall().second;
    return diameter;
}
};

template <typename Polytope, typename Point>
struct compute_diameter<BallIntersectPolytope<Polytope, Ball<Point>>>
{
template <typename NT>
static NT compute(BallIntersectPolytope<Polytope, Ball<Point>> const& P)
{
    NT diameter = NT(2) * P.radius();
    return diameter;
}
};

template <typename Point>
struct compute_diameter<ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>>
{
template <typename NT>
static NT compute(ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>> const& P)
{
    typedef typename ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>::VT VT;
    typedef typename ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>::MT MT;
    typedef HPolytope<Point> Hpolytope;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RandomNumberGenerator;
    PushBackWalkPolicy push_back_policy;
    typedef typename BCDHRWalk::template Walk
            <
                    Hpolytope,
                    RandomNumberGenerator
            > BCdhrWalk;
    typedef BoundaryRandomPointGenerator<BCdhrWalk> BCdhrRandomPointGenerator;

    MT G = P.get_T().transpose();
    MT AG = P.get_mat()*G;
    int k = G.cols(), d = P.dimension();
    MT eyes1(k, 2*k);
    eyes1 << MT::Identity(k,k), NT(-1) * MT::Identity(k,k);
    MT M1(k, 4*k);
    M1 << AG.transpose(), eyes1;
    MT M = M1.transpose();
    VT b = P.get_vec();

    VT bb(4*k);
    for (int i = 0; i < 4*k; ++i) bb(i) = (i < 2*k) ? b(i) : 1.0;

    Hpolytope HP;
    HP.init(d, M, bb);

    RandomNumberGenerator rng(HP.dimension());

    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall = HP.ComputeInnerBall();
    Point q = InnerBall.first;
    BCdhrRandomPointGenerator::apply(HP, q, 4*d*d, 1,
                                     randPoints, push_back_policy, rng);
    typename std::list<Point>::iterator rpit=randPoints.begin();
    NT max_norm = NT(0), iter_norm;
    for ( ; rpit!=randPoints.end(); rpit++) {
        iter_norm = (G*(*rpit).getCoefficients()).norm();
        if (iter_norm > max_norm) max_norm = iter_norm;
    }
    NT diameter = NT(2) * max_norm;
    return diameter;
}
};


// Billiard walk for uniform distribution

struct BilliardWalk
{
    BilliardWalk(double L)
            :   param(L, true)
    {}

    BilliardWalk()
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
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef HPolytope<Point> Hpolytope;
    typedef Zonotope<Point> zonotope;
    typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        initialize(P, p, rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      Point& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T = rng.sample_urdist() * _Len;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            T = rng.sample_urdist() * _Len;
            _v = GetDirection<Point>::apply(n, rng);
            Point p0 = _p;
            int it = 0;
            while (it < 50*n)
            {
                auto pbpair = P.line_positive_intersect(_p, _v, _lambdas,
                                                        _Av, _lambda_prev);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == 50*n){
                _p = p0;
            }
        }
        p = _p;
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }

private :

    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope const& P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _p = p;
        _v = GetDirection<Point>::apply(n, rng);

        NT T = rng.sample_urdist() * _Len;
        Point p0 = _p;
        int it = 0;

        std::pair<NT, int> pbpair
                = P.line_positive_intersect(_p, _v, _lambdas, _Av);
        if (T <= pbpair.first) {
            _p += (T * _v);
            _lambda_prev = T;
            return;
        }
        _lambda_prev = dl * pbpair.first;
        _p += (_lambda_prev * _v);
        T -= _lambda_prev;
        P.compute_reflection(_v, _p, pbpair.second);

        while (it <= 50*n)
        {
            std::pair<NT, int> pbpair
                    = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                break;
            }else if (it == 50*n) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                _p += (_lambda_prev * _v);
                break;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
        //if (it == 30*n) _p = p0;

    }

    NT _Len;
    Point _p;
    Point _v;
    NT _lambda_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};

};








#endif // RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP
