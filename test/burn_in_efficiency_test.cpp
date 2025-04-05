// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos
// Copyright (c) 2024 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <fstream>
#include <iostream>
#include <chrono>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "misc/misc.h"
#include "random_walks/random_walks.hpp"
#include "generators/known_polytope_generators.h"
#include "sampling/sampling.hpp"

// Test to compare the efficiency of the old and new burn-in implementations
TEST_CASE("burn_in_efficiency") {
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef HPolytope<Point> Polytope;
    
    // Create a simple polytope for testing
    Polytope P = generate_cube<Polytope>(3, 1.0, false);
    
    // Parameters for sampling
    unsigned int walkL = 10;
    unsigned int numpoints = 1000;
    unsigned int nburns = 1000;
    unsigned int d = P.dimension();
    
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;
    
    // Get a starting point
    P.get_chebyshev_center(StartingPoint);
    
    // Test with old implementation (using PushBackWalkPolicy and clear)
    auto start_old = std::chrono::high_resolution_clock::now();
    
    // Create a copy of the original sampling function with the old implementation
    typedef DikinWalk WalkType;
    typedef RandomPointGenerator<WalkType> RandomPointGenerator;
    PushBackWalkPolicy push_back_policy;
    
    // Old implementation
    if (nburns > 0) {
        RandomPointGenerator::apply(P, StartingPoint, nburns, walkL, randPoints,
                                    push_back_policy, rng);
        randPoints.clear();
    }
    RandomPointGenerator::apply(P, StartingPoint, numpoints, walkL, randPoints,
                                push_back_policy, rng);
    
    auto end_old = std::chrono::high_resolution_clock::now();
    auto duration_old = std::chrono::duration_cast<std::chrono::milliseconds>(end_old - start_old).count();
    
    // Clear for next test
    randPoints.clear();
    
    // Test with new implementation (using NoOpWalkPolicy)
    auto start_new = std::chrono::high_resolution_clock::now();
    
    // New implementation
    NoOpWalkPolicy no_op_policy;
    if (nburns > 0) {
        RandomPointGenerator::apply(P, StartingPoint, nburns, walkL, randPoints,
                                    no_op_policy, rng);
    }
    RandomPointGenerator::apply(P, StartingPoint, numpoints, walkL, randPoints,
                                push_back_policy, rng);
    
    auto end_new = std::chrono::high_resolution_clock::now();
    auto duration_new = std::chrono::duration_cast<std::chrono::milliseconds>(end_new - start_new).count();
    
    // Verify that both implementations produce the same number of points
    CHECK(randPoints.size() == numpoints);
    
    // Print the results
    std::cout << "Old implementation duration: " << duration_old << " ms" << std::endl;
    std::cout << "New implementation duration: " << duration_new << " ms" << std::endl;
    std::cout << "Efficiency improvement: " << (duration_old - duration_new) / (double)duration_old * 100 << "%" << std::endl;
    
    // The new implementation should be at least as fast as the old one
    CHECK(duration_new <= duration_old);
} 