// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLE_POINTS_HPP
#define SAMPLE_POINTS_HPP

struct UniformDistribution {};
struct GaussianDistribution {};

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
void sample_points(Polytope const& P,
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

        for (unsigned int i = 0; i < rnum; ++i)
        {
            walk.apply(P, p, walk_len, rng);
            samples.push_back(p);
        }
    }
    else if constexpr (std::is_same<WalkType, GaussianBallWalk>::value
               && std::is_same<Distribution, GaussianDistribution>::value)
    {
        //TODO
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