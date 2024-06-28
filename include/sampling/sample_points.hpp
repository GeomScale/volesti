// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLE_POINTS_HPP
#define SAMPLE_POINTS_HPP

#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_problem.h"
#include "sampling/sampling.hpp"

struct UniformDistribution
{
    using CartesianPoint = point<Cartesian<double>>;
    using Func = ZeroScalarFunctor<CartesianPoint>;
    using Grad = ZeroFunctor<CartesianPoint>;
    using Hess = ZeroFunctor<CartesianPoint>;
};

struct SphericalGaussianDistribution
{
    SphericalGaussianDistribution()
        :   variance(1.0)
    {}
    SphericalGaussianDistribution(double _variance)
        :   variance(_variance)
    {}
    double variance;
};

struct GaussianDistribution
{
    using CartesianEllipsoid = Ellipsoid<point<Cartesian<double>>>;

    GaussianDistribution() {}

    GaussianDistribution(CartesianEllipsoid _ellipsoid)
        :   ellipsoid(_ellipsoid)
    {}
    CartesianEllipsoid ellipsoid;
};

struct ExponentialDistribution
{
    using CartesianPoint = point<Cartesian<double>>;

    ExponentialDistribution() {}

    ExponentialDistribution(CartesianPoint _c, double _T)
        :   c(_c)
        ,   T(_T)
    {}
    CartesianPoint c;
    double T;
};


struct LogConcaveDistribution
{
    using CartesianPoint = point<Cartesian<double>>;

    template <typename GradientFunctor, typename FunctionFunctor>
    LogConcaveDistribution(GradientFunctor g, FunctionFunctor f, double _L)
        : grad(g)
        , func(f)
        , L(_L)
    {}

    template <typename GradientFunctor, typename FunctionFunctor, typename HessianFunctor>
    LogConcaveDistribution(GradientFunctor g, FunctionFunctor f, HessianFunctor h, double _L)
        : grad_point(g)
        , func(f)
        , hess(h)
        , L(_L)
    {}

    std::function<CartesianPoint(unsigned int const&, std::vector<CartesianPoint> const&, double const&)> grad;
    std::function<CartesianPoint(CartesianPoint const&)> grad_point;
    std::function<double(CartesianPoint const&)> func;
    std::function<CartesianPoint(CartesianPoint const&)> hess;
    double L;
};

namespace detail
{

template <typename T1, typename T2>
struct AlwaysFalse : std::false_type {};

template <typename WalkType, typename Distribution>
struct sample_points
{
    static_assert(AlwaysFalse<WalkType, Distribution>::value, "Not implemented for this combination of walk and distribution.");
};

}

template
<
    typename Polytope,
    typename Point,
    typename WalkType,
    typename Distribution,
    typename RandomNumberGenerator,
    typename PointList
>
void sample_points(Polytope& P, // TODO: make it a const&
                   Point const& starting_point,
                   WalkType const& walk_with_parameters,
                   Distribution const& distribution,
                   RandomNumberGenerator& rng,
                   unsigned int const& walk_len,
                   unsigned int const& rnum,
                   unsigned int const& nburns,
                   PointList& samples)
{
    if constexpr ((std::is_same<WalkType, AcceleratedBilliardWalk>::value
                || std::is_same<WalkType, BallWalk>::value
                || std::is_same<WalkType, BilliardWalk>::value
                || std::is_same<WalkType, CDHRWalk>::value
                || std::is_same<WalkType, DikinWalk>::value
                || std::is_same<WalkType, JohnWalk>::value
                || std::is_same<WalkType, RDHRWalk>::value
                || std::is_same<WalkType, VaidyaWalk>::value)
            && std::is_same<Distribution, UniformDistribution>::value)
    {
        typename WalkType::template Walk<Polytope, RandomNumberGenerator>
            walk(P, starting_point, rng, walk_with_parameters.param);

        Point p = starting_point;

        // accelerated billiard walk variants in some cases need to compute a good starting point
        // otherwise the walk stuck in the first given point (e.g. with a cube and its center of mass)
        if constexpr (std::is_same<WalkType, AcceleratedBilliardWalk>::value)
            walk.template get_starting_point(P, p, p, 10, rng);

        for (unsigned int i = 0; i < nburns; ++i)
        {
            walk.apply(P, p, walk_len, rng);
        }

        samples.reserve(rnum);
        for (unsigned int i = 0; i < rnum; ++i)
        {
            walk.apply(P, p, walk_len, rng);
            samples.push_back(p);
        }
    }
    else if constexpr ((std::is_same<WalkType, GaussianBallWalk>::value
                     || std::is_same<WalkType, GaussianCDHRWalk>::value
                     || std::is_same<WalkType, GaussianHamiltonianMonteCarloExactWalk>::value
                     || std::is_same<WalkType, GaussianRDHRWalk>::value)
                    && std::is_same<Distribution, SphericalGaussianDistribution>::value)
    {
        typename WalkType::template Walk<Polytope, RandomNumberGenerator>
            walk(P, starting_point, distribution.variance, rng, walk_with_parameters.param);

        Point p = starting_point;

        for (unsigned int i = 0; i < nburns; ++i)
        {
            walk.apply(P, p, distribution.variance, walk_len, rng);
        }

        samples.reserve(rnum);
        for (unsigned int i = 0; i < rnum; ++i)
        {
            walk.apply(P, p, distribution.variance, walk_len, rng);
            samples.push_back(p);
        }
    }
    else if constexpr (std::is_same<WalkType, GaussianAcceleratedBilliardWalk>::value
                    && std::is_same<Distribution, GaussianDistribution>::value)
    {
        typename WalkType::template Walk<Polytope, RandomNumberGenerator>
            walk(P, starting_point, distribution.ellipsoid, rng, walk_with_parameters.param);

        Point p = starting_point;

        for (unsigned int i = 0; i < nburns; ++i)
        {
            walk.apply(P, p, distribution.ellipsoid, walk_len, rng);
        }

        samples.reserve(rnum);
        for (unsigned int i = 0; i < rnum; ++i)
        {
            walk.apply(P, p, distribution.ellipsoid, walk_len, rng);
            samples.push_back(p);
        }
    }
    else if constexpr (std::is_same<WalkType, ExponentialHamiltonianMonteCarloExactWalk>::value
                    && std::is_same<Distribution, ExponentialDistribution>::value)
    {
        typename WalkType::template Walk<Polytope, RandomNumberGenerator>
            walk(P, starting_point, distribution.c, distribution.T, rng, walk_with_parameters.param);

        Point p = starting_point;

        for (unsigned int i = 0; i < nburns; ++i)
        {
            walk.apply(P, p, walk_len, rng);
        }

        samples.reserve(rnum);
        for (unsigned int i = 0; i < rnum; ++i)
        {
            walk.apply(P, p, walk_len, rng);
            samples.push_back(p);
        }
    }
    else if constexpr ((std::is_same<WalkType, HamiltonianMonteCarloWalk>::value
                     || std::is_same<WalkType, NutsHamiltonianMonteCarloWalk>::value
                     || std::is_same<WalkType, UnderdampedLangevinWalk>::value)
                    && std::is_same<Distribution, LogConcaveDistribution>::value)
    {
        using HPolytope = typename std::remove_const<Polytope>::type;

        using Solver = LeapfrogODESolver<Point, double, HPolytope, decltype(distribution.grad)>;

        std::vector<Point> xs;
        unsigned int i = 0;
        double t = 1.0;

        typename WalkType::parameters
        <
            double,
            decltype(distribution.grad)
        > hmc_params(distribution.L, P.dimension());

        Point walk_p = starting_point; //TODO: avoid the copy
        auto g = distribution.grad;
        auto f = distribution.func;

        HPolytope P_copy = P; //TODO: avoid the copy
        typename WalkType::template Walk
        <
            Point, HPolytope, RandomNumberGenerator,
            decltype(distribution.grad), decltype(distribution.func), Solver
        >
        walk(&P_copy, walk_p, g, f, hmc_params);

        Point p = starting_point;

        //TODO: burnin for nuts

        samples.reserve(rnum);
        for (int i = 0; i < rnum; i++)
        {
            walk.apply(rng, walk_len);
            samples.push_back(walk.x);
        }
    }
    else if constexpr ((std::is_same<WalkType, CRHMCWalk>::value)
                    && std::is_same<Distribution, LogConcaveDistribution>::value)
    {
        using HPolytope = typename std::remove_const<Polytope>::type;
        HPolytope HP = P; //TODO: avoid the copy

        constexpr int simdLen = 8; //TODO: input parameter
        using NT = double;

        int dimension = HP.dimension();

        auto G = distribution.grad_point;
        auto F = distribution.func;
        auto H = distribution.hess;

        using NegativeLogprobFunctor = decltype(distribution.func);
        using NegativeGradientFunctor = decltype(distribution.grad_point);
        using HessianFunctor = decltype(distribution.hess);

        using MatrixType = typename HPolytope::MT;
        using Input = crhmc_input
            <
                MatrixType,
                Point,
                NegativeLogprobFunctor,
                NegativeGradientFunctor,
                HessianFunctor
            >;
        Input input = convert2crhmc_input
            <
                Input, HPolytope, NegativeLogprobFunctor,
                NegativeGradientFunctor, HessianFunctor
            >(HP, F, G, H);

        using CrhmcProblem = crhmc_problem<Point, Input>;
        CrhmcProblem problem = CrhmcProblem(input);
        if(problem.terminate){return;}

        using Solver = ImplicitMidpointODESolver
            <
                Point, NT, CrhmcProblem,
                decltype(distribution.grad_point), simdLen
            >;

        using Walk = typename WalkType::template Walk
            <
                    Point,
                    CrhmcProblem,
                    RandomNumberGenerator,
                    NegativeGradientFunctor,
                    NegativeLogprobFunctor,
                    Solver
            >;
        using WalkParams = typename WalkType::template parameters
            <
                    NT,
                    NegativeGradientFunctor
            >;
        Point p = Point(problem.center);
        problem.options.simdLen=simdLen;
        WalkParams params(distribution.L, p.dimension(), problem.options);

        // TODO: pass a usef defined eta to bypass the one created by the walk using L
        //if (distribution.eta > 0) {
        //    params.eta = distribution.eta;
        // }

        Walk walk(problem, p, input.df, input.f, params);

        for (unsigned int i = 0; i < nburns; ++i)
        {
            walk.apply(rng, walk_len);
        }

        samples.reserve(rnum);
        bool raw_output = false; // TODO: check this
        for (unsigned int i = 0; i < std::ceil((float)rnum/simdLen); ++i)
        {
            walk.apply(rng, walk_len);
            if (walk.P.terminate) {return;}
            auto x = raw_output ? walk.x : walk.getPoints();

            if ((i + 1) * simdLen > rnum)
            {
              for (int j = 0; j < rnum-simdLen*i; j++)
              {
                samples.push_back(Point(x.col(j)));
              }
              break;
            }
            for (int j = 0; j < x.cols(); j++)
            {
              samples.push_back(Point(x.col(j)));
            }
        }
    }

    else
    {
        static_assert(detail::AlwaysFalse<WalkType, Distribution>::value,
            "Not implemented for this combination of walk and distribution.");
    }
}


// this is a wrapper function for storing the samples in an Eigen matrix
// it is inefficient in many ways (e.g. unneeded copies) and should be rewritten
template
<
    typename Polytope,
    typename Point,
    typename WalkType,
    typename Distribution,
    typename RandomNumberGenerator,
    typename NT
>
void sample_points(Polytope const& P,
                   Point const& starting_point,
                   WalkType const& walk_with_parameters,
                   Distribution const& distribution,
                   RandomNumberGenerator& rng,
                   unsigned int const& walk_len,
                   unsigned int const& rnum,
                   unsigned int const& nburns,
                   Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>& samples)
{
    std::vector<Point> random_point(1);
    auto p = starting_point;
    sample_points(P, p, walk_with_parameters, distribution, rng, walk_len, 0, nburns, random_point);
    random_point.clear();

    for (unsigned int i = 0; i < rnum; ++i)
    {
        sample_points(P, p, walk_with_parameters, distribution, rng, walk_len, 1, 0, random_point);
        auto point = random_point[0];
        samples.col(i) = point.getCoefficients();
        random_point.clear();
    }

}

#endif //SAMPLE_POINTS_HPP