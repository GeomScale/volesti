// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPE_GENERATORS_H
#define POLYTOPE_GENERATORS_H

#include <time.h>

template <class Polytope>
Polytope gen_cube(int dim, bool Vpoly) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    MT A;
    VT b;
    int m;

    if (!Vpoly) {

        A.resize(2 * dim, dim);
        b.resize(2 * dim);
        for (int i = 0; i < dim; ++i) {
            b(i) = 1.0;
            for (int j = 0; j < dim; ++j) {
                if (i == j) {
                    A(i, j) = 1.0;
                } else {
                    A(i, j) = 0.0;
                }
            }
        }
        for (int i = 0; i < dim; ++i) {
            b(i + dim) = 1.0;
            for (int j = 0; j < dim; ++j) {
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
        for(int i=0; i<m; ++i){
            b(i) = 1;
            int k=i, j=0;
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
    Polytope P;
    P.init(dim, A, b);

    return P;
}


template <class Polytope>
Polytope gen_cross(int dim, bool Vpoly) {

    int m;
    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    Polytope P;
    MT A;
    VT b;
    if (!Vpoly) {

        m = 2 << (dim - 1);
        A.resize(m, dim);
        b.resize(m);
        for(int i=0; i<m; ++i){
            b(i) = 1;
            int k=i, j=0;
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
    } else {
        A.resize(2 * dim, dim);
        b.resize(2 * dim);

        for(int i=0; i<dim; ++i){
            b(i) = 1.0;
            for(int j=0; j<dim; ++j){
                if(i==j) {
                    A(i,j) = 1.0;
                } else {
                    A(i,j) = 0.0;
                }
            }
        }
        for(int i=0; i<dim; ++i){
            b(i + dim) = 1.0;
            for(int j=0; j<dim; ++j){
                if(i==j) {
                    A(i + dim, j) = -1.0;
                } else {
                    A(i + dim, j) = 0.0;
                }
            }
        }
    }
    P.init(dim, A, b);
    return P;
}


template <class Polytope>
Polytope gen_simplex(int dim, bool Vpoly){
    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    A.resize(dim+1, dim);
    VT b;
    b.resize(dim+1);

    for(int i=0; i<dim; ++i){
        b(i) = 0;
        for(int j=0; j<dim; ++j){
            if(i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    b(dim) = 1;
    for(int j=0; j<dim; ++j){
        if (!Vpoly) {
            A(dim, j) = -1.0;
        } else {
            A(dim, j) = 0.0;
        }
    }
    Polytope P;
    P.init(dim, A, b);

    return P;
}


template <class Polytope>
Polytope gen_prod_simplex(int dim, bool Vpoly = false){

    if (Vpoly) {
        std::cout<<"Only prod simplices in H-representation can be generated.."<<std::endl;
        exit(-1);
    }

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    VT b;
    A.resize(2 * dim + 2, 2 * dim);
    b.resize(2 * dim + 2);
    Polytope P;

    //first simplex
    for(int i=0; i<dim; ++i){
        b(i) = 0.0;
        for(int j=0; j<dim; ++j) {
            A(i, j) = 0.0;
        }
        for(int j=0; j<dim; ++j){
            if(i==j) {
                A(i, j + dim) = 1.0;
            }else{
                A(i, j + dim) = 0.0;
            }
        }
    }

    b(dim) = 1.0;
    for(int j=0; j<dim; ++j) {
        A(dim, j) = 0.0;
    }
    for(int j=0; j<dim; ++j){
        A(dim, j + dim) = -1.0;
    }

    //second simplex
    for(int i=0; i<dim; ++i){
        b(dim + 1 + i) = 0.0;
        for(int j=0; j<dim; ++j){
            if(i==j) {
                A(dim + 1 + i, j) = 1.0;
            }else {
                A(dim + 1 + i, j) = 0.0;
            }
        }
        for(int j=0; j<dim; ++j) {
            A(dim + 1 + i, j + dim) = 0.0;
        }
    }
    b(2 * dim +1) = 1.0;
    for(int j=0; j<dim; ++j) {
        A(2 * dim +1, j) = -1.0;
    }
    for(int j=0; j<dim; ++j) {
        A(2 * dim +1, j + dim) = 0.0;
    }

    P.init(2 * dim, A, b);
    return P;
}


template <class Polytope>
Polytope gen_skinny_cube(int dim, bool Vpoly = false) {
    if (Vpoly) {
        std::cout<<"Only skinny cubes in H-representation can be generated.."<<std::endl;
        exit(-1);
    }

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    A.resize(2 * dim, dim);
    VT b;
    b.resize(2 * dim);

    for(int i=0; i<dim; ++i){
        if (i==0) {
            b(i) = 100.0;
        } else {
            b(i) = 1.0;
        }
        for(int j=0; j<dim; ++j){
            if(i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    for(int i=0; i<dim; ++i){
        if (i==0) {
            b(i + dim) = 100.0;
        } else {
            b(i + dim) = 1.0;
        }
        for(int j=0; j<dim; ++j){
            if(i==j) {
                A(i + dim, j) = -1.0;
            } else {
                A(i + dim, j) = 0.0;
            }
        }
    }
    Polytope P;
    P.init(dim, A, b);

    return P;
}


template <class Polytope>
Polytope gen_prod_simplex(int dim, int m){

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;

    MT A;
    VT b;
    A.resize(dim, m);
    b.resize(dim);
    Polytope P;

    srand (time(NULL));

    for(int i=1; i<m; ++i){
        b(dim) = 1000.0;
        for(int j=0; j<dim; ++j){
            if(std::rand()%2==1) {
                A(i,j) = -std::rand()%1000;
            } else {
                A(i,j) = std::rand()%1000;
            }
        }
    }
    P.init(dim, A, b);
    return P;
}


template <class Polytope, class RNGType>
Polytope gen_zonotope(int dim, int m) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);

    NT sum;
    MT A;
    VT b;
    A.resize(m, dim);
    b.resize(m);
    Polytope P;

    for (int i = 0; i < dim; ++i) {
        b(i) = 0.0;
        for (int j = 0; j < dim; ++j) {
            if (i==j) {
                A(i,j) = 1.0;
            } else {
                A(i,j) = 0.0;
            }
        }
    }

    for (int i = dim; i < m; ++i) {
        b(i) = 0.0;
        sum = 0.0;
        for (int j = 0; j < dim; ++j) {
            A(i,j) = rdist(rng);
            sum += A(i,j) * A(i,j);
        }
        sum = 1.0 / std::sqrt(sum);
        for (int j = 0; j < dim; ++j) {
            A(i,j) = A(i,j) * sum;
        }
    }

    P.init(dim, A, b);
    return P;
}


#endif
