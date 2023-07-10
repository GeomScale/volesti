// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef V_POLYTOPES_GEN_H
#define V_POLYTOPES_GEN_H

#include <exception>

#ifndef isnan
  using std::isnan;
#endif

template <class MT>
void removeRow(MT &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

/// Generates a random V-polytope
/// @tparam Polytope polytope type
/// @tparam RNGType RNGType type
template <class Polytope, class RNGType>
Polytope random_vpoly(unsigned int dim, unsigned int k, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PointType PointType;
    typedef PointType Point;

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }
    boost::normal_distribution<> rdist(0,1);

    typename std::vector<NT>::iterator pit;
    MT V(k, dim);
    unsigned int j;


    std::vector<NT> Xs(dim,0);
    NT normal = NT(0);

    for (unsigned int i = 0; i < k; ++i) {

        normal = NT(0);
        for (unsigned int i=0; i<dim; i++) {
            Xs[i] = rdist(rng);
            normal += Xs[i] * Xs[i];
        }
        normal = 1.0 / std::sqrt(normal);

        for (unsigned int i=0; i<dim; i++) {
            Xs[i] = Xs[i] * normal;
        }

        for (unsigned int j=0; j<dim; j++) {
            V(i,j) = Xs[j];
        }
    }


    VT b = VT::Ones(k);
    return Polytope(dim, V, b);
}

/// Generates a random V-polytope inside a cube
/// @tparam Polytope polytope type
/// @tparam RNGType RNGType type
template <class Polytope, class RNGType>
Polytope random_vpoly_incube(unsigned int d, unsigned int k, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PointType PointType;
    typedef PointType Point;

    REAL *conv_mem;
    int *colno_mem;

    conv_mem = (REAL *) malloc(k * sizeof(*conv_mem));
    colno_mem = (int *) malloc(k * sizeof(*colno_mem));

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    Point p(d);
    typename std::vector<NT>::iterator pit;
    MT V(k, d);
    unsigned int j, count_row,it=0;
    std::vector<int> indices;

    VT b = VT::Ones(k);

    for (unsigned int i = 0; i < k; ++i) {
        for (int j = 0; j < d; ++j) {
            V(i, j) = urdist1(rng);
        }
    }
    if(k==d+1){
        return Polytope(d, V, b);
    }

    MT V2(k,d);
    V2 = V;
    indices.clear();
    while(it<20) {
        V.resize(V2.rows(), d);
        V = V2;
        for (int i = 0; i < indices.size(); ++i) {
            V.conservativeResize(V.rows()+1, d);
            for (int j = 0; j < d; ++j) {
                V(V.rows()-1, j) = urdist1(rng);
            }
        }
        indices.clear();
        V2.resize(k, d);
        V2 = V;

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < d; ++j) {
                p.set_coord(j, V(i, j));
            }
            removeRow(V2, i);
            if (memLP_Vpoly(V2, p, conv_mem, colno_mem)){
                indices.push_back(i);
            }
            V2.resize(k, d);
            V2 = V;
        }
        if (indices.size()==0) {
            return Polytope(d, V, b);
        }
        V2.resize(k - indices.size(), d);
        count_row =0;
        for (int i = 0; i < k; ++i) {
            if(std::find(indices.begin(), indices.end(), i) != indices.end()) {
                continue;
            } else {
                for (int j = 0; j < d; ++j) V2(count_row, j) = V(i,j);
                count_row++;
            }
        }
        it++;
    }


    free(colno_mem);
    free(conv_mem);

    return Polytope(d, V2, VT::Ones(V2.rows()));
//    return VP;

}

#endif
