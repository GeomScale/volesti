// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2022 Vissarion Fisikopoulos
// Copyright (c) 2018-2022 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_COMPUTE_DIAMETER_HPP
#define RANDOM_WALKS_COMPUTE_DIAMETER_HPP

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#ifndef DISABLE_LPSOLVE
    #include "convex_bodies/vpolytope.h"
    #include "convex_bodies/vpolyintersectvpoly.h"
    #include "convex_bodies/zpolytope.h"
    #include "convex_bodies/zonoIntersecthpoly.h"
#endif
#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/ellipsoid.h"


template <typename GenericPolytope>
struct compute_diameter
{
    template <typename NT>
    static NT compute(GenericPolytope) {return NT(0);}
};


template <typename Point>
struct compute_diameter<HPolytope<Point>>
{
    template <typename NT>
    static NT compute(HPolytope<Point> &P)
    {
        return NT(2) * std::sqrt(NT(P.dimension())) * P.InnerBall().second;
    }
};

template <typename Point>
struct compute_diameter<Spectrahedron<Point>>
{
    template <typename NT>
    static NT compute(Spectrahedron<Point> &P)
    {
        std::pair<Point, NT> inner_ball = P.ComputeInnerBall();
        return NT(6) * NT(P.dimension()) * inner_ball.second;
    }
};

template <typename Point>
struct compute_diameter<CorrelationSpectrahedron<Point>>
{
    template <typename NT>
    static NT compute(CorrelationSpectrahedron<Point> &P)
    {
        std::pair<Point, NT> inner_ball = P.getInnerBall();
        return NT(P.dimension()) * inner_ball.second;
    }
};

template <typename Point>
struct compute_diameter<CorrelationSpectrahedron_MT<Point>>
{
    template <typename NT>
    static NT compute(CorrelationSpectrahedron_MT<Point> &P)
    {
        std::pair<Point, NT> inner_ball = P.getInnerBall();
        return NT(P.dimension()) * inner_ball.second;
    }
};

#ifndef DISABLE_LPSOLVE
template <typename Point>
struct compute_diameter<VPolytope<Point>>
{
    template <typename NT>
    static NT compute(VPolytope<Point> &P)
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
    static NT compute(Zonotope<Point> &P)
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

        return NT(2) * (V.transpose() * x0).norm();
    }
};

template <typename Point, typename RandomNumberGenerator>
struct compute_diameter<IntersectionOfVpoly<VPolytope<Point>, RandomNumberGenerator>>
{
    template <typename NT>
    static NT compute(IntersectionOfVpoly<VPolytope<Point>, RandomNumberGenerator> &P)
    {
        return NT(2) * NT(P.dimension()) * P.InnerBall().second;
    }
};

template <typename Point>
struct compute_diameter<ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>>>
{
    template <typename NT>
    static NT compute(ZonoIntersectHPoly<Zonotope<Point>, HPolytope<Point>> &P)
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

        Hpolytope HP(d, M, bb);

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
        return NT(2) * max_norm;
    }
};

#endif

template <typename Polytope, typename Point>
struct compute_diameter<BallIntersectPolytope<Polytope, Ball<Point>>>
{
    template <typename NT>
    static NT compute(BallIntersectPolytope<Polytope, Ball<Point>> &P)
    {
        return NT(2) * P.radius();
    }
};

template <typename Point>
struct compute_diameter<OrderPolytope<Point>>
{
    template <typename NT>
    static NT compute(OrderPolytope<Point>& P)
    {
        return std::sqrt(NT(P.dimension()));
    }
};

template <typename Point>
struct compute_diameter<BallIntersectPolytope<OrderPolytope<Point>, Ellipsoid<Point> > >
{
    template <typename NT>
    static NT compute(BallIntersectPolytope<OrderPolytope<Point>, Ellipsoid<Point>>& P)
    {
        NT polytope_diameter = std::sqrt(NT(P.dimension()));
        return std::min(polytope_diameter, (NT(2) * P.radius()));
    }
};

#endif // RANDOM_WALKS_COMPUTE_DIAMETER_HPP
