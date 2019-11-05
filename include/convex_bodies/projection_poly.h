// VolEsti (volume computation and sampling library)

// Copyright (c) 2019 Vissarion Fisikopoulos
// Copyright (c) 2019 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef PROJECTION_POLY_H
#define PROJECTION_POLY_H

#include <limits>

#include <iostream>
#include "solve_lp.h"
#include "projection_oracles.h"

//min and max values for the Hit and Run functions

// V-Polytope class
template <class Point, class  RNGType>
class ProjPoly{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef RNGType rngtype;

private:
    MT T;  //matrix V. Each row contains a vertex
    VT b;  // vector b that contains first column of ine file
    MT A;
    unsigned int _d;  //dimension
    REAL *conv_comb, *row;
    int *colno;
    NT maxNT = std::numeric_limits<NT>::max();

public:
    ProjPoly() {}

    // return dimension
    unsigned int dimension() {
        return _d;
    }


    // this function returns 0. The main sampler requests this function to set the length of lambdas vector
    int num_of_hyperplanes() {
        return 0;
    }


    // compute the number of facets of the cyclic polytope in dimension _d with the same number of vertices
    // this is an upper bound for the number of the facets from McMullen's Upper Bound Theorem
    unsigned int upper_bound_of_hyperplanes() {

        return 2*_d;
    }


    // return the number of vertices
    int num_of_vertices() {
        return T.rows();
    }


    // return the matrix V
    MT get_mat() {
        return T;
    }


    // return the vector b
    VT get_vec() {
        return b;
    }


    // change the matrix V
    void set_mat(MT V2) {
        T = V2;
    }


    // change the vector b
    void set_vec(VT b2) {
        b = b2;
    }


    // get a specific coeff of matrix V
    NT get_mat_coeff(unsigned int i, unsigned int j) {
        return T(i,j);
    }


    // get a specific coeff of vector b
    NT get_vec_coeff(unsigned int i) {
        return b(i);
    }


    // set a specific coeff of matrix V
    void put_mat_coeff(unsigned int i, unsigned int j, NT value) {
        T(i,j) = value;
    }


    // set a specific coeff of vector b
    void put_vec_coeff(unsigned int i, NT value) {
        b(i) = value;
    }

    void comp_diam(NT &diam) {
        return;
    }


    void init(unsigned int dim, MT _T, MT _A, VT _b) {
        _d = dim;
        T = _T;
        A = _A;
        b = _b;
        conv_comb = (REAL *) malloc((T.cols()+1) * sizeof(*conv_comb));
        colno = (int *) malloc((T.cols()+1) * sizeof(*colno));
        row = (REAL *) malloc((T.cols()+1) * sizeof(*row));
    }


    // Construct matrix V which contains the vertices row-wise
    void init(std::vector<std::vector<NT> > Pin) {
        _d = Pin[0][1] - 1;
        T.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                T(i - 1, j - 1) = Pin[i][j];
            }
        }
        conv_comb = (REAL *) malloc(Pin.size() * sizeof(*conv_comb));
        colno = (int *) malloc((T.rows()+1) * sizeof(*colno));
        row = (REAL *) malloc((T.rows()+1) * sizeof(*row));
    }


    // print polytope in input format
    void print() {
#ifdef VOLESTI_DEBUG
        std::cout << " " << V.rows() << " " << _d << " float" << std::endl;
#endif
        for (unsigned int i = 0; i < T.rows(); i++) {
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



    // pick d+1 random vertices until they define a full dimensional simplex and then
    // compute the chebychev ball of that simplex
    std::pair<Point,NT> ComputeInnerBall() {

        std::pair<Point, NT> che_up = ComputeChebychevBall<NT, Point>(A, b);
        b = b - A*Eigen::Map<VT>(&(che_up.first).get_coeffs()[0], (che_up.first).dimension());

        std::vector<NT> temp(_d,0);
        NT radius =  maxNT, min_plus;
        std::pair<NT,NT> min_max;
        Point center(_d);

        for (unsigned int i = 0; i < _d; ++i) {
            temp.assign(_d,0);
            temp[i] = 1.0;
            Point v(_d,temp.begin(), temp.end());
            min_max = intersect_double_line_proj_poly(T, A, b, center, v, conv_comb, row, colno);
            
            if (radius > min_max.first) radius = min_max.first;
            if (radius > -min_max.second) radius = -min_max.second;
        }

        radius = radius / std::sqrt(NT(_d));
        return std::pair<Point, NT> (center, radius);
    }


    // check if point p belongs to the convex hull of V-Polytope P
    int is_in(Point p) {
        if(memLP_proj_poly(T, A, b, p)){
            return -1;
        }
        return 0;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v) {

        return intersect_double_line_proj_poly(T, A, b, r, v, conv_comb, row, colno);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        return intersect_double_line_proj_poly(T, A, b, r, v, conv_comb, row, colno);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        return intersect_double_line_proj_poly(T, A, b, r, v, conv_comb, row, colno);
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {
        return std::pair<NT, int> (intersect_line_proj_poly(T, A, b, r, v, conv_comb, row, colno), 1);
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av,
                                               NT &lambda_prev) {
        return line_positive_intersect(r, v, Ar, Av);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());
        return intersect_double_line_proj_poly(T, A, b, r, v, conv_comb, row, colno);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());
        return intersect_double_line_proj_poly(T, A, b, r, v, conv_comb, row, colno);
    }


    // shift polytope by a point c
    void shift(VT c) {
       return;
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(MT _T) {
        T = _T*T;
    }

    Point get_mean_of_vertices() {
        return Point(_d);
    }

    NT get_max_vert_norm() {
        return 0.0;
    }

    // consider an upper bound for the number of facets of a V-polytope
    // for each facet consider a lower bound for the distance from the origin
    // useful for CG algorithm to get the first gaussian
    std::vector<NT> get_dists(NT radius) {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    void normalize() {}


    // in number_of_vertices<=20*dimension use the vertices for the rounding
    // otherwise you have to sample from the V-polytope
    template <class PointList>
    bool get_points_for_rounding (PointList &randPoints) {
        return false;
    }

    MT get_T() {
        return T;
    }

    void compute_reflection(Point &v, Point &p, int facet) {

        //Eigen::FullPivLU<MatrixXd> lu(A);
        //Eigen::MatrixXd A_null_space = lu.kernel();
        return;
    }

    void free_them_all() {
        free(row);
        free(colno);
        free(conv_comb);
    }

};

#endif
