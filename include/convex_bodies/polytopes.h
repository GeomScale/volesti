// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <iostream>

//min and max values for the Hit and Run functions
const NT maxNT = 1.79769e+308;
const NT minNT = -1.79769e+308;

// my H-polytope class
template <typename FT>
class Polytope{
private:
    Eigen::MatrixXd A; //matrix A
    Eigen::VectorXd b; // vector b, s.t.: Ax<=b
    int            _d; //dimension

public:
    Polytope() {}

    // constructor: cube(d)
    Polytope(int d): _d(d) {
        A.resize(2 * d, d);
        b.resize(2 * d);
        for (int i = 0; i < d; ++i) {
            b(i) = 1;
            for (int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i, j) = 1;
                } else {
                    A(i, j) = 0;
                }
            }
        }
        for (int i = 0; i < d; ++i) {
            b(i + d) = 1;
            for (int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i + d, j) = -1;
                } else {
                    A(i + d, j) = 0;
                }
            }
        }
    }

    // return dimension
    int dimension() {
        return _d;
    }

    // return the number of facets
    int num_of_hyperplanes() {
        return A.rows();
    }

    // return the matrix A
    Eigen::MatrixXd get_eigen_mat() {
        return A;
    }

    // return the vector b
    Eigen::VectorXd get_eigen_vec() {
        return b;
    }

    // chenge the matrix A
    void set_eigen_mat(Eigen::MatrixXd A2) {
        A = A2;
    }

    // change the vector b
    void set_eigen_vec(Eigen::VectorXd b2) {
        b = b2;
    }

    // get a single coefficient of the matrix or the vector
    FT get_coeff(int i, int j) {
        if (j == 0) {
            return b(i);
        }
        return A(i, j - 1);
    }

    // change a single coefficient of the matrix or the vector
    void put_coeff(int i, int j, FT value) {
        if (j == 0) {
            b(i) = value;
        } else {
            A(i, j - 1) = value;
        }
    }

    // default initialize: cube(d)
    void init(int d) {
        A.resize(2 * d, d);
        b.resize(2 * d);
        for (int i = 0; i < d; ++i) {
            b(i) = 1;
            for (int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i, j) = 1;
                } else {
                    A(i, j) = 0;
                }
            }
        }
        for (int i = 0; i < d; ++i) {
            b(i + d) = 1;
            for (int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i + d, j) = -1;
                } else {
                    A(i + d, j) = 0;
                }
            }
        }
    }

    //define matrix A and vector b, s.t. Ax<=b and the dimension
    void init(std::vector<std::vector<FT> > Pin) {
        _d = Pin[0][1] - 1;
        A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (int j = 1; j < _d + 1; j++) {
                A(i - 1, j - 1) = Pin[i][j];
            }
        }
    }

    // print polytope in input format
    void print() {
        std::cout << " " << A.rows() << " " << _d + 1 << " float" << std::endl;
        for (int i = 0; i < A.rows(); i++) {
            for (int j = 0; j < _d; j++) {
                std::cout << A(i, j) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
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
    int is_in(Point p) {
        FT sum;
        int m = A.rows(), i, j;
        for (i = 0; i < m; i++) {
            sum = b(i);
            for (j = 0; j < _d; j++) {
                sum -= A(i, j) * p[j];
            }
            if (sum < FT(0)) { //Check if corresponding hyperplane is violated
                return 0;
            }
        }
        return -1;
    }


    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point,FT> chebyshev_center() {

        std::pair <Point,FT> res;
        res = solveLP(A, b, _d);  //lpSolve lib for the linear program
        return res;
    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<FT,FT> line_intersect(Point r,
                                          Point v) {

        FT lamda = 0, min_plus = FT(maxNT), max_minus = FT(minNT);
        FT sum_nom, sum_denom;
        int i, j, m = num_of_hyperplanes();
        typename std::vector<FT>::iterator rit, vit;

        for (i = 0; i < m; i++) {
            sum_nom = b(i);
            sum_denom = FT(0);
            j = 0;
            rit = r.iter_begin();
            vit = v.iter_begin();
            for ( ; rit != r.iter_end(); rit++, vit++, j++){
                sum_nom -= A(i, j) * (*rit);
                sum_denom += A(i, j) * (*vit);
            }
            if (sum_denom == FT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = sum_nom / sum_denom;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }
        return std::pair<FT, FT>(min_plus, max_minus);
    }


    //First coordinate ray intersecting convex polytope
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          int rand_coord,
                                          std::vector<FT> &lamdas) {

        FT lamda = 0, min_plus = FT(maxNT), max_minus = FT(minNT);
        FT sum_nom, sum_denom;
        int i, j, m = num_of_hyperplanes();
        typename std::vector<FT>::iterator rit;

        for (i = 0; i < m; i++) {
            sum_nom = b(i);
            sum_denom = A(i, rand_coord);
            rit = r.iter_begin();
            j = 0;
            for (; rit != r.iter_end(); rit++, j++) {
                sum_nom -= A(i, j) * (*rit);
            }
            lamdas[i] = sum_nom;
            if (sum_denom == FT(0)) {
                //std::cout<<"div0"<<sum_denom<<std::endl;
                ;
            } else {
                lamda = sum_nom * (1 / sum_denom);
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
        }
        return std::pair<FT, FT>(min_plus, max_minus);
    }


    //Not the first coordinate ray intersecting convex
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<FT> &lamdas) {

        typename std::vector<FT>::iterator lamdait = lamdas.begin();
        FT lamda = 0, min_plus = FT(maxNT), max_minus = FT(minNT);
        FT sum_nom, sum_denom, c_rand_coord, c_rand_coord_prev;
        int i, j, m = num_of_hyperplanes();
        typename std::vector<FT>::iterator rit;

        for (i = 0; i < m; i++) {
            sum_denom = b(i);
            c_rand_coord = A(i, rand_coord);
            c_rand_coord_prev = A(i, rand_coord_prev);

            *lamdait = *lamdait
                       + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
            if (c_rand_coord == FT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = (*lamdait) / c_rand_coord;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            ++lamdait;
        }
        return std::pair<FT, FT>(min_plus, max_minus);
    }


    //Apply linear transformation, of square matrix T, in H-polytope P:= Ax<=b
    void linear_transformIt(Eigen::MatrixXd T) {
        A = A * T;
    }



};

#endif
