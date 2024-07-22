// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ORDER_POLYTOPES_GEN_H
#define ORDER_POLYTOPES_GEN_H

#include <chrono>
#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "misc/misc.h"
#include "misc/poset.h"

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "generators/boost_random_number_generator.hpp"

#include "convex_bodies/orderpolytope.h"
#include "convex_bodies/hpolytope.h"


// Instances taken from: https://github.com/ttalvitie/le-counting-practice
static const std::unordered_map<std::string, std::string> instances =
{
    {"bipartite_0.5_008_0", R"(0 0 0 0 1 0 1 0
                               0 0 0 0 1 0 0 0
                               0 0 0 0 1 1 0 1
                               0 0 0 0 1 0 1 0
                               0 0 0 0 0 0 0 0
                               0 0 0 0 0 0 0 0
                               0 0 0 0 0 0 0 0
                               0 0 0 0 0 0 0 0)"},

    {"bayesiannetwork_andes_008_0", R"(0 0 0 0 0 0 0 0
                                       1 0 0 0 0 0 0 0
                                       0 1 0 0 0 0 0 0
                                       0 0 0 0 0 0 0 0
                                       0 0 1 1 0 0 0 0
                                       0 0 0 0 1 0 0 0
                                       0 0 0 0 0 0 0 0
                                       0 0 0 0 1 1 1 0)"},

};

// generates a Polytope from a poset
/// @tparam Polytope Type of returned polytope
template <class Polytope>
Polytope get_orderpoly(Poset const &poset) {
    typedef typename Polytope::PointType Point;

    OrderPolytope<Point> OP(poset);
    if constexpr (std::is_same< Polytope, OrderPolytope<Point> >::value ) {
        return OP;
    } else if constexpr (std::is_same<Polytope, HPolytope<Point> >::value ){
        Polytope HP(OP.dimension(), OP.get_dense_mat(), OP.get_vec());
        return HP;
    } else {
        throw "Unable to generate an Order Polytope of requested type";
    }
}

// generates an Order Polytope from an instance name
// Instances taken from: https://github.com/ttalvitie/le-counting-practice
/// @tparam Polytope Type of returned polytope
template <class Polytope>
Polytope generate_orderpoly(std::string& instance_name) {
    std::stringstream in_ss(instances.at(instance_name));
    Poset poset = read_poset_from_file_adj_matrix(in_ss).second;
    return get_orderpoly<Polytope>(poset);
}

// Generates a cube as an Order Polytope
/// @tparam Polytope Type of returned polytope
template <class Polytope>
Polytope generate_cube_orderpoly(unsigned int dim) {
    typedef typename Poset::RV RV;

    RV order_relations;
    Poset poset(dim, order_relations);
    return get_orderpoly<Polytope>(poset);
}

// Generates a random Order Polytope with given dimension and number of facets
/// @tparam Polytope Type of returned polytope
/// @tparam RNGType RNGType Type
template <class Polytope, typename NT>
Polytope random_orderpoly(unsigned int dim, unsigned int m, int seed = std::numeric_limits<int>::signaling_NaN()) {

    typedef typename Poset::RV RV;

    int rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    if (!isnan(seed)) {
        rng_seed = seed;
    }

    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNG;
    RNG rng(dim);
    rng.set_seed(rng_seed);


    std::vector<unsigned int> order(dim);
    for(int i = 0; i < dim; ++i) {
        order[i] = i;
    }
    boost::mt19937 shuffle_rng(rng_seed);
    std::shuffle(order.begin(), order.end(), shuffle_rng);


    RV order_relations;
    for(int i = 0; i < m - 2 * dim; ++i) {
        unsigned int x = rng.sample_uidist();
        unsigned int y = rng.sample_uidist();
        while(x == y) {
            y = rng.sample_uidist();
        }
        if(x > y)
            std::swap(x, y);
        order_relations.push_back(std::make_pair(order[x], order[y]));
    }


    Poset poset(dim, order_relations);
    return get_orderpoly<Polytope>(poset);
}

#endif
