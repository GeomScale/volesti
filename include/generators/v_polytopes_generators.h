// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef V_POLYTOPES_GEN_H
#define V_POLYTOPES_GEN_H

#include <exception>

template <class MT>
void removeRow(MT &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

template <class Polytope, class RNGType>
Polytope random_vpoly(unsigned int dim, unsigned int k, double seed = std::numeric_limits<double>::signaling_NaN()) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PointType PointType;
    typedef PointType Point;

    unsigned rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(rng_seed);
    if (!std::isnan(seed)) {
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

    Polytope VP;
    VT b = VT::Ones(k);
    VP.init(dim, V, b);

    return VP;

}


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
    if (!std::isnan(seed)) {
        unsigned rng_seed = seed;
        rng.seed(rng_seed);
    }
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    Point p(d);
    typename std::vector<NT>::iterator pit;
    MT V(k, d);
    unsigned int j, count_row,it=0;
    std::vector<int> indices;
    Polytope VP;
    VT b = VT::Ones(k);

    for (unsigned int i = 0; i < k; ++i) {
        for (int j = 0; j < d; ++j) {
            V(i, j) = urdist1(rng);
        }
    }
    if(k==d+1){
        VP.init(d, V, b);
        return VP;
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
            VP.init(d, V, b);
            return VP;
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

    VP.init(d, V2, VT::Ones(V2.rows()));
    free(colno_mem);
    free(conv_mem);

    return VP;

}

#endif
