// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef HPOLYTOPE_H
#define HPOLYTOPE_H

#include <limits>

#include <iostream>
#include "solve_lp.h"

//min and max values for the Hit and Run functions


// H-polytope class
template <typename Point>
class HPolytope{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    //using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    //typedef RowMatrixXd MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    MT A; //matrix A
    VT b; // vector b, s.t.: Ax<=b
    unsigned int            _d; //dimension
    //NT maxNT = 1.79769e+308;
    //NT minNT = -1.79769e+308;
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();

public:
    HPolytope() {}

    // constructor: cube(d)
    HPolytope(unsigned int d): _d(d) {
        A.resize(2 * d, d);
        b.resize(2 * d);
        for (unsigned int i = 0; i < d; ++i) {
            b(i) = 1;
            for (unsigned int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i, j) = 1;
                } else {
                    A(i, j) = 0;
                }
            }
        }
        for (unsigned int i = 0; i < d; ++i) {
            b(i + d) = 1;
            for (unsigned int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i + d, j) = -1;
                } else {
                    A(i + d, j) = 0;
                }
            }
        }
    }


    // return dimension
    unsigned int dimension() const {
        return _d;
    }


    // return the number of facets
    int num_of_hyperplanes() const {
        return A.rows();
    }


    // return the matrix A
    MT get_mat() const {
        return A;
    }


    // return the vector b
    VT get_vec() const {
        return b;
    }


    // change the matrix A
    void set_mat(const MT &A2) {
        A = A2;
    }


    // change the vector b
    void set_vec(const VT &b2) {
        b = b2;
    }


    // set a specific coeff of matrix A
    NT get_mat_coeff(const unsigned int &i, const unsigned int &j) const {
        return A(i,j);
    }


    // get a spesific coeff of vector b
    NT get_vec_coeff(const unsigned int &i) const {
        return b(i);
    }


    // get a specific coeff of matrix A
    void put_mat_coeff(const unsigned int &i, const unsigned int &j, const NT &value) {
        A(i,j) = value;
    }


    // set a spesific coeff of vector b
    void put_vec_coeff(const unsigned int &i, const NT &value) {
        b(i) = value;
    }


    Point get_mean_of_vertices() const {
        return Point(_d);
    }


    NT get_max_vert_norm() const {
        return 0.0;
    }

    void comp_diam(NT &diam, const NT &cheb_rad) {
        if(cheb_rad < 0.0) {
            diam = 2.0 * std::sqrt(NT(_d)) * ComputeInnerBall().second;
        } else {
            diam = 2.0 * std::sqrt(NT(_d)) * cheb_rad;
        }
    }

    void init(const unsigned int dim, const MT &_A, const VT &_b) {
        _d = dim;
        A = _A;
        b = _b;
    }

    //define matrix A and vector b, s.t. Ax<=b and the dimension
    void init(const std::vector<std::vector<NT> > &Pin) {
        _d = Pin[0][1] - 1;
        A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                A(i - 1, j - 1) = -Pin[i][j];
            }
        }
    }


    // print polytope in input format
    void print() {
#ifdef VOLESTI_DEBUG
        std::cout << " " << A.rows() << " " << _d + 1 << " float" << std::endl;
#endif
        for (unsigned int i = 0; i < A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                #ifdef VOLESTI_DEBUG
                std::cout << -A(i, j) << " ";
                #endif
            }
            #ifdef VOLESTI_DEBUG
            std::cout << "<= " << b(i) << std::endl;
            #endif
        }
    }


    // Compute the reduced row echelon form
    // used to transofm {Ax=b,x>=0} to {A'x'<=b'}
    // e.g. Birkhoff polytopes
    /*
    // Todo: change the implementation in order to use eigen matrix and vector.
    int rref(){
        to_reduced_row_echelon_form(_A);
        std::vector<int> zeros(_d+1,0);
        std::vector<int> ones(_d+1,0);
        std::vector<int> zerorow(_A.size(),0);
        for (int i = 0; i < _A.size(); ++i)
        {
            for (int j = 0; j < _d+1; ++j){
                if ( _A[i][j] == double(0)){
                    ++zeros[j];
                    ++zerorow[i];
                }
                if ( _A[i][j] == double(1)){
                    ++ones[j];
                }
            }
        }
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            int j =0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ){
                if(zeros[j]==_A.size()-1 && ones[j]==1)
                    (*mit).erase(lit);
                else{ //reverse sign in all but the first column
                    if(lit!=mit->end()-1) *lit = (-1)*(*lit);
                    ++lit;
                }
                ++j;
            }
        }
        //swap last and first columns
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            double temp=*(mit->begin());
            *(mit->begin())=*(mit->end()-1);
            *(mit->end()-1)=temp;
        }
        //delete zero rows
        for (typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ){
            int zero=0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit){
                if(*lit==double(0)) ++zero;
            }
            if(zero==(*mit).size())
                _A.erase(mit);
            else
                ++mit;
        }
        //update _d
        _d=(_A[0]).size();
        // add unit vectors
        for(int i=1;i<_d;++i){
            std::vector<double> e(_d,0);
            e[i]=1;
            _A.push_back(e);
        }
        // _d should equals the dimension
        _d=_d-1;
        return 1;
    }*/

    
    //Check if Point p is in H-polytope P:= Ax<=b
    int is_in(const Point &p) const {
        NT sum;
        int m = A.rows();
        for (int i = 0; i < m; i++) {
            sum = b(i) - A.row(i) * p.getCoefficients();

            //Check if corresponding hyperplane is violated
            if (sum < NT(0)) return 0;
        }
        return -1;
    }


    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point,NT> ComputeInnerBall() {

        //lpSolve lib for the linear program
        return ComputeChebychevBall<NT, Point>(A, b);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point &r, Point &v) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes();
        viterator rit, vit;

        for (int i = 0; i < m; i++) {
            sum_nom = b(i) - A.row(i) * r.getCoefficients();
            sum_denom = A.row(i) * v.getCoefficients();

            if (sum_denom == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = sum_nom / sum_denom;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point &r, Point &v, std::vector<NT> &Ar,
            std::vector<NT> &Av, bool pos = false) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom, mult;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;
        viterator rit, vit, Ariter = Ar.begin(), Aviter = Av.begin();

        for (int i = 0; i < m; i++, ++Ariter, ++Aviter) {
            sum_nom = -A.row(i) * r.getCoefficients();
            sum_denom = A.row(i) * v.getCoefficients();

            (*Ariter) = -sum_nom;
            (*Aviter) = sum_denom;
            sum_nom += b(i);
            if (sum_denom == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = sum_nom / sum_denom;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, std::vector<NT> &Ar,
            std::vector<NT> &Av, const NT &lambda_prev, bool pos = false) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom, mult;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;
        viterator vit, Ariter = Ar.begin(), Aviter = Av.begin();

        for (int i = 0; i < m; i++, ++Ariter, ++Aviter) {
            (*Ariter) += lambda_prev * (*Aviter);
            sum_nom = b(i) - (*Ariter);
            sum_denom = A.row(i) * v.getCoefficients();

            (*Aviter) = sum_denom;
            if (sum_denom == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = sum_nom / sum_denom;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos) facet = i;
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point &r, Point &v, std::vector<NT> &Ar,
            std::vector<NT> &Av) {
        return line_intersect(r, v, Ar, Av, true);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point &r, Point &v, std::vector<NT> &Ar,
            std::vector<NT> &Av, const NT &lambda_prev) {
        return line_intersect(r, v, Ar, Av, lambda_prev, true);
    }


    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point &r, const unsigned int &rand_coord,
                                          std::vector<NT> &lamdas) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom;
        unsigned int j;
        int m = num_of_hyperplanes();
        viterator rit;

        for (int i = 0; i < m; i++) {
            sum_nom = b(i) - A.row(i)*r.getCoefficients();
            sum_denom = A(i, rand_coord);
            lamdas[i] = sum_nom;
            if (sum_denom == NT(0)) {
                //std::cout<<"div0"<<sum_denom<<std::endl;
                ;
            } else {
                lamda = sum_nom * (1 / sum_denom);
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const Point &r_prev,
                                          const unsigned int rand_coord,
                                          const unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        viterator lamdait = lamdas.begin();
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom, c_rand_coord, c_rand_coord_prev;
        int m = num_of_hyperplanes();

        for (int i = 0; i < m; i++) {
            sum_denom = b(i);
            c_rand_coord = A(i, rand_coord);
            c_rand_coord_prev = A(i, rand_coord_prev);

            *lamdait = *lamdait + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
            if (c_rand_coord == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = (*lamdait) / c_rand_coord;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            ++lamdait;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    void linear_transformIt(const MT &T) {
        A = A * T;
    }


    // shift polytope by a point c
    void shift(const VT &c){
        b = b - A*c;
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(const NT &radius){
        unsigned int i=0;
        std::vector <NT> dists(num_of_hyperplanes(), NT(0));
        typename std::vector<NT>::iterator disit = dists.begin();
        for ( ; disit!=dists.end(); disit++, i++)
            *disit = b(i) / A.row(i).norm();

        return dists;
    }


    // no points given for the rounding, you have to sample from the polytope
    template <typename T>
    bool get_points_for_rounding (const T &randPoints) {
        return false;
    }

    MT get_T() const {
        return A;
    }

    void normalize() {

        NT row_norm;
        for (int i = 0; i < num_of_hyperplanes(); ++i) {
            row_norm = A.row(i).norm();
            A.row(i) = A.row(i) / row_norm;
            b(i) = b(i) / row_norm;
        }

    }

    void compute_reflection(Point &v, const Point &p, const int facet) {

        VT a = A.row(facet);
        Point s(_d, std::vector<NT>(&a[0], a.data()+a.cols()*a.rows()));
        s = ((-2.0 * v.dot(s)) * s);
        v = s + v;

    }

    void free_them_all() {}

};

#endif
