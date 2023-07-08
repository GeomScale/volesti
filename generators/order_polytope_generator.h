// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ORDER_POLYTOPES_GEN_H
#define ORDER_POLYTOPES_GEN_H

#include <sstream>
#include <unordered_map>
#include "misc.h"
#include "misc/poset.h"


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

// generates an Order Polytope from an instance name
// Instances taken from: https://github.com/ttalvitie/le-counting-practice
/// @tparam Polytope Type of returned polytope
template <class Polytope>
Polytope generate_orderpoly(std::string& instance_name) {
    std::stringstream in_ss(instances.at(instance_name));
    Poset poset = read_poset_from_file_adj_matrix(in_ss).second;
    return Polytope(poset);
}

// Generates a cube as an Order Polytope
/// @tparam Polytope Type of returned polytope
template <class Polytope>
Polytope generate_cube_orderpoly(unsigned int dim) {
    typedef typename Poset::RV RV;

    RV order_relations;
    Poset poset(dim, order_relations);
    Polytope OP(poset);
    return OP;
}

#endif
