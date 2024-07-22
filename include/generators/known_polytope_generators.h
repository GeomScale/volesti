// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef KNOWN_POLYTOPE_GENERATORS_H
#define KNOWN_POLYTOPE_GENERATORS_H

#include <exception>

#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"

/// This function generates a hypercube of given dimension
/// The result can be either in V-representation (Vpoly=true) or in H-representation (V-poly-false)
/// @tparam Polytope Type of returned polytope
template <class Polytope>
Polytope generate_cube(const unsigned int& dim, const bool& Vpoly,typename Polytope::NT  scale=1) {

    typedef typename Eigen::Matrix<typename Polytope::NT, Eigen::Dynamic, Eigen::Dynamic>    MT;
    typedef typename Polytope::VT    VT;
    MT A;
    VT b;
    unsigned int m;

    if (!Vpoly) {

        A.resize(2 * dim, dim);
        b.resize(2 * dim);
        for (unsigned int i = 0; i < dim; ++i) {
            b(i) = scale;
            for (unsigned int j = 0; j < dim; ++j) {
                if (i == j) {
                    A(i, j) = 1.0;
                } else {
                    A(i, j) = 0.0;
                }
            }
        }
        for (unsigned int i = 0; i < dim; ++i) {
            b(i + dim) = scale;
            for (unsigned int j = 0; j < dim; ++j) {
                if (i == j) {
                    A(i + dim, j) = -1.0;
                } else {
                    A(i + dim, j) = 0.0;
                }
            }
        }
    } else {

        m = 2 << (dim - 1);
        A.resize(m, dim);
        b.resize(m);
        for(unsigned int i=0; i<m; ++i){
            b(i) = 1;
            unsigned int k=i, j=0;
            while(k!=0){
                int bit = k % 2;
                if(bit == 0) {
                    A(i,j) = -1.0;
                }else {
                    A(i,j) = 1.0;
                }
                k = k >> 1;
                ++j;
            }
            for(; j<dim; ++j) {
                A(i,j) = -1.0;
            }
        }
    }

    return Polytope(dim, A, b);
}


/// This function generates a crosspolytope of given dimension
/// The result can be either in V-representation (Vpoly=true) or in H-representation (V-poly-false)
/// @tparam Polytope Type of returned polytope
template <typename Polytope>
Polytope generate_cross(const unsigned int &dim, const bool &Vpoly) {

    unsigned int m;
    typedef typename Eigen::Matrix<typename Polytope::NT, Eigen::Dynamic, Eigen::Dynamic>    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    VT b;
    if (!Vpoly) {

        m = 2 << (dim - 1);
        A.resize(m, dim);
        b.resize(m);
        for(unsigned int i=0; i<m; ++i){
            b(i) = 1;
            unsigned int k=i, j=0;
            while(k!=0){
                unsigned int bit = k % 2;
                if(bit == 0) {
                    A(i,j) = -1.0;
                }else {
                    A(i,j) = 1.0;
                }
                k = k >> 1;
                ++j;
            }
            for(; j<dim; ++j) {
                A(i,j) = -1.0;
            }
        }
    } else {
        A.resize(2 * dim, dim);
        b.resize(2 * dim);

        for(unsigned int i=0; i<dim; ++i){
            b(i) = 1.0;
            for(unsigned int j=0; j<dim; ++j){
                if(i==j) {
                    A(i,j) = 1.0;
                } else {
                    A(i,j) = 0.0;
                }
            }
        }
        for(unsigned int i=0; i<dim; ++i){
            b(i + dim) = 1.0;
            for(unsigned int j=0; j<dim; ++j){
                if(i==j) {
                    A(i + dim, j) = -1.0;
                } else {
                    A(i + dim, j) = 0.0;
                }
            }
        }
    }
    return Polytope(dim, A, b);
}


/// This function generates a simplex of given dimension
/// The result can be either in V-representation (Vpoly=true) or in H-representation (V-poly-false)
/// @tparam Polytope Type of returned polytope
template <typename Polytope>
Polytope generate_simplex(const unsigned int &dim, const bool &Vpoly){
    typedef typename Eigen::Matrix<typename Polytope::NT, Eigen::Dynamic, Eigen::Dynamic>    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    A.resize(dim+1, dim);
    VT b;
    b.resize(dim+1);

    for(unsigned int i=0; i<dim; ++i){
        if (!Vpoly) {
            b(i) = 0;
        } else {
            b(i) = 1;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    b(dim) = 1;
    for(unsigned int j=0; j<dim; ++j){
        if (!Vpoly) {
            A(dim, j) = -1.0;
        } else {
            A(dim, j) = 0.0;
        }
    }

    return Polytope(dim, A, b);
}


/// This function generates a product of simplices of given dimension
/// The result can be either in V-representation (Vpoly=true) or in H-representation (V-poly-false)
/// @tparam Polytope Type of returned polytope
template <typename Polytope>
Polytope generate_prod_simplex(const unsigned int &dim, bool Vpoly = false){

    Polytope Perr;
    try
    {
        if(Vpoly) throw false;
    }
    catch (bool e) {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Only prod simplices in H-representation can be generated.."<<std::endl;
        #endif
        return Perr;
    }

    typedef typename Eigen::Matrix<typename Polytope::NT, Eigen::Dynamic, Eigen::Dynamic>    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    VT b;
    A.resize(2 * dim + 2, 2 * dim);
    b.resize(2 * dim + 2);


    //first simplex
    for(unsigned int i=0; i<dim; ++i){
        b(i) = 0.0;
        for(unsigned int j=0; j<dim; ++j) {
            A(i, j) = 0.0;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i, j + dim) = 1.0;
            }else{
                A(i, j + dim) = 0.0;
            }
        }
    }

    b(dim) = 1.0;
    for(unsigned int j=0; j<dim; ++j) {
        A(dim, j) = 0.0;
    }
    for(unsigned int j=0; j<dim; ++j){
        A(dim, j + dim) = -1.0;
    }

    //second simplex
    for(unsigned int i=0; i<dim; ++i){
        b(dim + 1 + i) = 0.0;
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(dim + 1 + i, j) = 1.0;
            }else {
                A(dim + 1 + i, j) = 0.0;
            }
        }
        for(unsigned int j=0; j<dim; ++j) {
            A(dim + 1 + i, j + dim) = 0.0;
        }
    }
    b(2 * dim +1) = 1.0;
    for(unsigned int j=0; j<dim; ++j) {
        A(2 * dim +1, j) = -1.0;
    }
    for(unsigned int j=0; j<dim; ++j) {
        A(2 * dim +1, j + dim) = 0.0;
    }

    return Polytope(2 * dim, A, b);
}


/// This function generates a skinny cube of given dimension
/// The result can be either in V-representation (Vpoly=true) or in H-representation (V-poly-false)
/// @tparam Polytope Type of returned polytope
template <typename Polytope>
Polytope generate_skinny_cube(const unsigned int &dim, bool Vpoly = false) {

    Polytope Perr;
    try
    {
        if(Vpoly) throw false;
    }
    catch (bool e) {
        #ifdef VOLESTI_DEBUG
        std::cout<<"Only prod simplices in H-representation can be generated.."<<std::endl;
        #endif
        return Perr;
    }

    typedef typename Eigen::Matrix<typename Polytope::NT, Eigen::Dynamic, Eigen::Dynamic>    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    A.resize(2 * dim, dim);
    VT b;
    b.resize(2 * dim);

    for(unsigned int i=0; i<dim; ++i){
        if (i==0) {
            b(i) = 100.0;
        } else {
            b(i) = 1.0;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    for(unsigned int i=0; i<dim; ++i){
        if (i==0) {
            b(i + dim) = 100.0;
        } else {
            b(i + dim) = 1.0;
        }
        for(unsigned int j=0; j<dim; ++j){
            if(i==j) {
                A(i + dim, j) = -1.0;
            } else {
                A(i + dim, j) = 0.0;
            }
        }
    }
    return Polytope(dim, A, b);
}

/// This function generates the Birkhoff polytope of given type n
/// The Birkhoff polytope also called the assignment polytope or the polytope of doubly stochastic matrices.
/// @tparam Polytope Type of returned polytope
template <typename Polytope>
Polytope generate_birkhoff(unsigned int const& n) {

    unsigned int m = n * n;
    unsigned int d = n * n - 2 * n + 1;

    typedef typename Eigen::Matrix<typename Polytope::NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef typename Polytope::VT VT;

    MT A = MT::Zero(m, d);
    VT b(m);

    b(d) = -1.0 * int(n - 2);

    for (int i = 0; i < d; ++i) {
        A(d, i) = -1;
    }

    for (int i = 0; i < d; ++i) {
        b(i) = 0;
        A(i, i) = -1;
    }

    for (int i = d+1; i < d+1+n-1; ++i) {
        b(i) = 1;
        for (int counter = 0; counter < n-1; ++counter) {
            A(i, counter * (n-1) + (i-d-1)) = 1;
        }
    }

    for (int i = d+n; i < m; ++i) {
        b(i) = 1;
        for (int counter = 0; counter < n-1; ++counter) {
            A(i, counter + (i-d-n) * (n-1)) = 1;
        }
    }

    Polytope P(d, A, b);

    return P;
}

#endif
