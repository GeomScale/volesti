// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

// Count the number of linear extensions of a poset through volume
// approximation of the corresponding order polytope.
//
// Key identity (Stanley, 1986):
//   #LinearExtensions(P) = n! * Vol(OrderPolytope(P))
//
// This module provides a clean API for computing this quantity
// using any of volesti's volume estimation algorithms.

#ifndef COUNT_LINEAR_EXTENSIONS_HPP
#define COUNT_LINEAR_EXTENSIONS_HPP

#include <cmath>
#include "misc/poset.h"
#include "convex_bodies/orderpolytope.h"

#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"


// Compute n! (factorial) as a double
template <typename NT>
NT factorial(unsigned int n)
{
    NT result = NT(1);
    for (unsigned int i = 2; i <= n; ++i)
        result *= NT(i);
    return result;
}


///
/// Count linear extensions of a poset using volume estimation of the order polytope.
///
/// @tparam WalkTypePolicy   The walk type policy (e.g., CDHRWalk, AcceleratedBilliardWalk, etc.)
/// @tparam RandomNumberGenerator  The RNG type
/// @tparam NT               Numeric type (defaults to double)
///
/// @param poset         The partial order set
/// @param error         Upper bound for volume approximation error (default: 0.1)
/// @param walk_length   Walk length per sample (default: 1)
///
/// @return Approximate number of linear extensions
///
template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename NT = double
>
NT count_linear_extensions(Poset const& poset,
                           NT const& error = 0.1,
                           unsigned int const& walk_length = 1)
{
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef OrderPolytope<Point> OrderPoly;

    // Build the order polytope from the poset
    OrderPoly OP(poset);
    unsigned int d = OP.dimension();

    RandomNumberGenerator rng(d);

    // Compute volume of the order polytope using cooling balls (default)
    NT volume = volume_cooling_balls<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length).second;

    // #LE = n! * Vol(OrderPolytope)
    return factorial<NT>(d) * volume;
}


///
/// Count linear extensions using a specific volume algorithm.
///
/// @param algo  Algorithm choice: "SOB", "CG", or "CB" (default)
///
template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename NT = double
>
NT count_linear_extensions(Poset const& poset,
                           std::string const& algo,
                           NT const& error = 0.1,
                           unsigned int const& walk_length = 1)
{
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef OrderPolytope<Point> OrderPoly;

    OrderPoly OP(poset);
    unsigned int d = OP.dimension();

    RandomNumberGenerator rng(d);

    NT volume;
    if (algo == "SOB" || algo == "sob") {
        volume = volume_sequence_of_balls<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length);
    } else if (algo == "CG" || algo == "cg") {
        volume = volume_cooling_gaussians<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length);
    } else {
        // Default: Cooling Balls
        volume = volume_cooling_balls<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length).second;
    }

    return factorial<NT>(d) * volume;
}


///
/// Count linear extensions from a pre-built OrderPolytope.
/// Avoids re-constructing the polytope when called multiple times.
///
template
<
    typename WalkTypePolicy = CDHRWalk,
    typename RandomNumberGenerator = BoostRandomNumberGenerator<boost::mt19937, double>,
    typename Point,
    typename NT = double
>
NT count_linear_extensions_from_polytope(OrderPolytope<Point> & OP,
                                         std::string const& algo = "CB",
                                         NT const& error = 0.1,
                                         unsigned int const& walk_length = 1)
{
    unsigned int d = OP.dimension();
    RandomNumberGenerator rng(d);

    NT volume;
    if (algo == "SOB" || algo == "sob") {
        volume = volume_sequence_of_balls<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length);
    } else if (algo == "CG" || algo == "cg") {
        volume = volume_cooling_gaussians<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length);
    } else {
        volume = volume_cooling_balls<WalkTypePolicy, RandomNumberGenerator>(OP, rng, error, walk_length).second;
    }

    return factorial<NT>(d) * volume;
}


#endif // COUNT_LINEAR_EXTENSIONS_HPP
