// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ZPOLYTOPE_H
#define ZPOLYTOPE_H

#include <limits>

#include <iostream>
#include "solve_lp.h"

//min and max values for the Hit and Run functions

template <class Point>
class Zonotope {
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

private:
    MT V;  //matrix V. Each row contains a vertex
    VT b;  // vector b that contains first column of ine file
    unsigned int _d;  //dimension
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();
    REAL *conv_comb;
    MT sigma;
    MT Q0;
    MT T;

public:

    Zonotope() {}

    // return the dimension
    unsigned int dimension() {
        return _d;
    }


    // this function returns 0. The main sampler requests this function to set the length of lambdas vector
    int num_of_hyperplanes() {
        return 0;
    }


    // return the number of parallelopipeds. Used in get_dists fnction.
    unsigned int upper_bound_of_hyperplanes() {
        unsigned int m = V.rows(), d = _d, res;
        long double nom = 1.0, denom = 1.0, num_of_hyp = 0.0;

        for (unsigned int i = d+1 ; i <= m; ++i) {
            nom *= i;
        }
        for (unsigned int i = 1 ; i <= m-d; ++i) {
            denom *= i;
        }

        num_of_hyp = nom / denom;

        res = 2*_d;
        return res;
    }

    void compute_eigenvectors(MT G, bool norm1, bool norm2) {

        int k = G.cols();
        //std::cout<<"number of generators = "<<k<<std::endl;

        //MT ps = G.completeOrthogonalDecomposition().pseudoInverse();
        MT ps = G;//*G.transpose();
        sigma.resize(k,k);
        sigma = ps.transpose()*ps;
        //std::cout<<sigma<<std::endl;
        if (norm1) {
            sigma = (sigma + sigma.transpose()) / 2;
        }
        //std::cout<<"\n"<<sigma<<std::endl;
        //std::cout<<"\n";
        if (norm2) {
            for (int i1 = 0; i1 < k; ++i1) {
                sigma(i1,i1) = sigma(i1,i1) + 0.00000001;
            }
        }
        //MT sigma2 = sigma.inverse();
        Eigen::SelfAdjointEigenSolver<MT> es(sigma);

        MT D = es.eigenvalues().asDiagonal();
        MT Q2 = es.eigenvectors();

        //std::cout<<Q2<<"\n"<<std::endl;
        //std::cout<<es.eigenvalues()<<"\n"<<std::endl;
        //std::cout<<D<<"\n"<<std::endl;

        Q0.resize(k,k-_d);
        int count=0;
        for (int i = 0; i < k; ++i) {
            if(es.eigenvalues()[i]<0.0000001) {
                for (int j = 0; j < k; ++j) {
                    Q0(j, count) = Q2(j, i);
                }
                count++;
            }
        }
        Eigen::JacobiSVD<MT> svd(Q0, Eigen::ComputeFullU | Eigen::ComputeFullV);
        MT T2 = svd.matrixU().transpose();
        T.resize(_d,k);
        for (int i = k-_d; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                T(i-k+_d,j) = T2(i,j);
            }
        }
        //std::cout<<G<<"\n"<<std::endl;
        //std::cout<<Q0<<"\n"<<std::endl;
        for (int i1 = 0; i1 < k; ++i1) {
            sigma(i1,i1) = sigma(i1,i1) + 0.00000001;
        }

    }

    MT get_T() {
        return T;
    }

    MT get_Q0(){
        return Q0;
    }

    MT get_sigma() {
        return sigma;
    }

    // return the number of vertices
    int num_of_vertices() {
        return V.rows();
    }


    // return the number of generators
    int num_of_generators() {
        return V.rows();
    }


    // return the matrix V
    MT get_mat() {
        return V;
    }


    // return the vector b
    VT get_vec() {
        return b;
    }


    // change the matrix V
    void set_mat(MT V2) {
        V = V2;
    }


    // change the vector b
    void set_vec(VT b2) {
        b = b2;
    }


    // get a specific coeff of matrix V
    NT get_mat_coeff(unsigned int i, unsigned int j) {
        return V(i,j);
    }


    // get a specific coeff of vector b
    NT get_vec_coeff(unsigned int i) {
        return b(i);
    }


    // set a specific coeff of matrix V
    void put_mat_coeff(unsigned int i, unsigned int j, NT value) {
        V(i,j) = value;
    }


    // set a specific coeff of vector b
    void put_vec_coeff(unsigned int i, NT value) {
        b(i) = value;
    }


    // define zonotope using Eigen matrix V. Vector b is neded in order the code to compatible with Hpolytope class
    void init(unsigned int dim, MT _V, VT _b) {
        _d = dim;
        V = _V;
        b = _b;
        conv_comb = (REAL *) malloc((V.rows()+1) * sizeof(*conv_comb));
        bool normalization1=true;
        bool normalization2=false;
        compute_eigenvectors(V.transpose(),normalization1,normalization2);
    }


    // Construct matrix V which contains the vertices row-wise
    void init(std::vector<std::vector<NT> > Pin) {
        _d = Pin[0][1] - 1;
        V.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                V(i - 1, j - 1) = Pin[i][j];
            }
        }
        conv_comb = (REAL *) malloc(Pin.size() * sizeof(*conv_comb));
        bool normalization1=true;
        bool normalization2=false;
        compute_eigenvectors(V.transpose(),normalization1,normalization2);
    }


    // print polytope in input format
    void print() {
#ifdef VOLESTI_DEBUG
        std::cout << " " << V.rows() << " " << _d << " float" << std::endl;
#endif
        for (unsigned int i = 0; i < V.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
#ifdef VOLESTI_DEBUG
                std::cout << V(i, j) << " ";
#endif
            }
#ifdef VOLESTI_DEBUG
            std::cout<<"\n";
#endif
        }
    }


    // check if point p belongs to the convex hull of V-Polytope P
    int is_in(Point p) {
        if(memLP_Zonotope(V, p)){
            return -1;
        }
        return 0;
    }


    // Compute an inner ball of the zonotope
    std::pair<Point,NT> ComputeInnerBall() {
        std::vector<NT> temp(_d,0);
        NT radius =  maxNT, min_plus;
        Point center(_d);

        for (unsigned int i = 0; i < _d; ++i) {
            temp.assign(_d,0);
            temp[i] = 1.0;
            Point v(_d,temp.begin(), temp.end());
            min_plus = intersect_line_Vpoly<NT>(V, center, v, conv_comb, false, true);
            if (min_plus < radius) radius = min_plus;
        }

        radius = radius / std::sqrt(NT(_d));
        return std::pair<Point, NT> (center, radius);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point r, Point v) {

        return intersect_line_zono<NT>(V, r, v);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        return intersect_line_zono<NT>(V, r, v);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        return intersect_line_zono<NT>(V, r, v);
    }

    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        std::pair<NT, int> vppair;
        vppair.first = intersect_line_Vpoly(V, r, v, conv_comb, false, true);
        vppair.second = 1;
        return vppair;
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av,
                                               NT &lambda_prev) {
        return line_positive_intersect(r, v, Ar, Av);
    }


    // Compute the intersection of a coordinate ray
    // with the Zonotope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d,0);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        return intersect_line_zono<NT>(V, r, v);

    }


    // Compute the intersection of a coordinate ray
    // with the Zonotope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d,0);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        return intersect_line_zono<NT>(V, r, v);
    }


    // shift polytope by a point c
    // vector c has to be always the zero vector
    void shift(VT c) {
        return;
    }


    // get number of parallelopipeds
    // for each parallelopiped consider a lower bound for the distance from the origin
    // useful for CV algorithm to get the first gaussian
    std::vector<NT> get_dists(NT radius) {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    // apply linear transformation, of square matrix T, to thr Zonotope
    void linear_transformIt(MT T) {
        MT V2 = T.inverse() * V.transpose();
        V = V2.transpose();
    }

    // return false to the rounding function
    // no points are given so they have o be sampled
    template <class T>
    bool get_points_for_rounding (T &randPoints) {
        return false;
    }

    void compute_reflection(Point &v, Point &p, int facet) {

        int count = 0, outvert;
        MT Fmat(_d-1,_d);
        NT e = 0.0000000001;
        for (int j = 0; j < num_of_generators(); ++j) {
            if (((1.0 - *(conv_comb + j) ) > e || (1.0 - *(conv_comb + j) ) > e*std::abs(*(conv_comb + j))) && ((1.0 + *(conv_comb + j) ) > e || (1.0 + *(conv_comb + j) ) > e*std::abs(*(conv_comb + j)))) {
                Fmat.row(count) = V.row(j);
                count++;
            } else {
                outvert = j;
            }
        }

        VT a = Fmat.fullPivLu().kernel();
        if (a.dot(V.row(outvert)) > 1.0) a = -a;
        a = a/a.norm();

        Point s(_d);
        for (int i = 0; i < _d; ++i) {
            s.set_coord(i, a(i));
        }
        s = ((-2.0 * v.dot(s)) * s);
        v = s + v;
    }

};

#endif
