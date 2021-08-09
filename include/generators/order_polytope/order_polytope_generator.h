#include <iostream>
#include "misc.h"
#include "misc/poset.h"

// generates an Order Polytope from an instance file containing the adjacency matrix
// of the corresponding DAG
template <class Polytope>
Polytope generate_orderpoly(std::string& filename) {
    Poset poset = read_poset_from_file_adj_matrix(filename).second;
    return Polytope(poset);
}

// generates a cube as an Order Polytope
template <class Polytope>
Polytope generate_cube_orderpoly(unsigned int dim) {
    typedef typename Poset::RV RV;

    RV order_relations;
    Poset poset(dim, order_relations);
    Polytope OP(poset);
    return OP;
}