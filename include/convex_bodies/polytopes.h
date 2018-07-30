// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <iostream>


// my H-polytope class
template <typename K>
class Polytope{
private:
    Eigen::MatrixXd A; //matrix A
    Eigen::VectorXd b; // vector b, s.t.: Ax<=b
    int            _d; //dimension

public:
    typedef K                    FT;
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

    int dimension() {
        return _d;
    }

    int num_of_hyperplanes() {
        return A.rows();
    }

    Eigen::MatrixXd get_eigen_mat() {
        return A;
    }

    Eigen::VectorXd get_eigen_vec() {
        return b;
    }

    int set_eigen_mat(Eigen::MatrixXd A2) {
        A = A2;
        return -1;
    }

    int set_eigen_vec(Eigen::VectorXd b2) {
        b = b2;
        return -1;
    }

    K get_coeff(int i, int j) {
        if (j == 0) {
            return b(i);
        }
        return A(i, j - 1);
    }

    void put_coeff(int i, int j, K value) {
        if (j == 0) {
            b(i) = value;
        } else {
            A(i, j - 1) = value;
        }
    }

    // default initialize: cube(d)
    int init(int d) {
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
        return 0;
    }

    int init(std::vector<std::vector<K> > Pin) {
        _d = Pin[0][1] - 1;
        //define matrix A and vector b, s.t. Ax<=b
        A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (int j = 1; j < _d + 1; j++) {
                A(i - 1, j - 1) = Pin[i][j];
            }
        }
        return 0;
    }

    // print polytope in input format
    int print() {
        std::cout << " " << A.rows() << " " << _d + 1 << " float" << std::endl;
        for (int i = 0; i < A.rows(); i++) {
            for (int j = 0; j < _d; j++) {
                std::cout << A(i, j) << " ";
            }
            std::cout << "<= " << b(i) << std::endl;
        }
        return 0;
    }


    // Compute the reduced row echelon form
    // used to transofm {Ax=b,x>=0} to {A'x'<=b'}
    // e.g. Birkhoff polytopes
    /*
    // Todo change the implementation in order to use eigen matrix.
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

    

    int is_in(Point p) {
        K sum;
        int m = A.rows(), i, j;
        for (i = 0; i < m; i++) {
            sum = b(i);
            for (j = 0; j < _d; j++) {
                sum -= A(i, j) * p[j];
            }
            if (sum < K(0)) {
                return 0;
            }
        }
        return -1;
    }


    std::pair<Point,double> chebyshev_center() {

        std::pair<Point, double> res;
        res = solveLP(A, b, _d);
        return res;

    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point r,
                                          Point v) {
        K lamda = 0;
        K min_plus = 0, max_minus = 0;
        K sum_nom, sum_denom;
        bool min_plus_not_set = true;
        bool max_minus_not_set = true;
        int m = num_of_hyperplanes();

        for (int i = 0; i < m; i++) {
            sum_nom = b(i);
            sum_denom = K(0);
            for (int j = 0; j < _d; j++) {
                sum_nom -= A(i, j) * r[j];
                sum_denom += A(i, j) * v[j];
            }
            if (sum_denom == K(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = sum_nom / sum_denom;
                if (min_plus_not_set && lamda > 0) {
                    min_plus = lamda;
                    min_plus_not_set = false;
                }
                if (max_minus_not_set && lamda < 0) {
                    max_minus = lamda;
                    max_minus_not_set = false;
                }
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          int rand_coord) {
            K lamda = 0;
            K min_plus = 0, max_minus = 0;
            bool min_plus_not_set = true;
            bool max_minus_not_set = true;
            K sum_nom, sum_denom;
            int m = num_of_hyperplanes();

            for (int i = 0; i < m; i++) {
                sum_nom = b(i);
                sum_denom = A(i, rand_coord);
                for (int j = 0; j < _d; j++) {
                    sum_nom -= A(i, j) * r[j];
                }
                if (sum_denom == K(0)) {
                    //std::cout<<"div0"<<sum_denom<<std::endl;
                    ;
                } else {
                    lamda = sum_nom * (1 / sum_denom);

                    if (min_plus_not_set && lamda > 0) {
                        min_plus = lamda;
                        min_plus_not_set = false;
                    }
                    if (max_minus_not_set && lamda < 0) {
                        max_minus = lamda;
                        max_minus_not_set = false;
                    }
                    if (lamda < min_plus && lamda > 0) min_plus = lamda;
                    if (lamda > max_minus && lamda < 0) max_minus = lamda;
                }
            }
            return std::pair<NT, NT>(min_plus, max_minus);
        }


    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init) {
            K lamda = 0;
            std::vector<NT>::iterator lamdait = lamdas.begin();

            K min_plus = 0, max_minus = 0;
            K sum_nom, sum_denom, c_rand_coord, c_rand_coord_prev;
            bool min_plus_not_set = true;
            bool max_minus_not_set = true;
            int mini, maxi, m = num_of_hyperplanes();

            if (init) { //first time compute the innerprod cit*rit
                for (int i = 0; i < m; i++) {
                    sum_nom = b(i);
                    sum_denom = A(i, rand_coord);
                    for (int j = 0; j < _d; j++) {
                        sum_nom -= A(i, j) * r[j];
                    }
                    lamdas[i] = sum_nom;
                    if (sum_denom == K(0)) {
                        //std::cout<<"div0"<<sum_denom<<std::endl;
                        ;
                    } else {
                        lamda = sum_nom * (1 / sum_denom);

                        if (min_plus_not_set && lamda > 0) {
                            min_plus = lamda;
                            min_plus_not_set = false;
                        }
                        if (max_minus_not_set && lamda < 0) {
                            max_minus = lamda;
                            max_minus_not_set = false;
                        }
                        if (lamda < min_plus && lamda > 0) min_plus = lamda;
                        if (lamda > max_minus && lamda < 0) max_minus = lamda;

                    }
                }
            } else {//only a few opers no innerprod
                for (int i = 0; i < m; i++) {
                    sum_denom = b(i);
                    c_rand_coord = A(i, rand_coord);
                    c_rand_coord_prev = A(i, rand_coord_prev);

                    *lamdait = *lamdait
                               + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
                    if (c_rand_coord == K(0)) {
                        //std::cout<<"div0"<<std::endl;
                        ;
                    } else {
                        lamda = (*lamdait) / c_rand_coord;

                        if (min_plus_not_set && lamda > 0) {
                            min_plus = lamda;
                            min_plus_not_set = false;
                        }
                        if (max_minus_not_set && lamda < 0) {
                            max_minus = lamda;
                            max_minus_not_set = false;
                        }
                        if (lamda < min_plus && lamda > 0) min_plus = lamda;
                        if (lamda > max_minus && lamda < 0) max_minus = lamda;

                    }
                    ++lamdait;
                }
            }
            return std::pair<NT, NT>(min_plus, max_minus);
        }

    
    int linear_transformIt(Eigen::MatrixXd T) {
        A = A * T;
        return 1;
    }



};

#endif
