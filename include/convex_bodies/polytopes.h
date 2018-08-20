// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <iostream>
#include "solve_lp.h"

//min and max values for the Hit and Run functions
const NT maxNT = 1.79769e+308;
const NT minNT = -1.79769e+308;

// H-polytope class
template <typename FT>
class HPolytope{
public:
    typedef Eigen::Matrix<FT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<FT,Eigen::Dynamic,1> VT;

private:
    MT A; //matrix A
    VT b; // vector b, s.t.: Ax<=b
    int            _d; //dimension
    typedef typename std::vector<FT>::iterator viterator;

public:
    HPolytope() {}

    // constructor: cube(d)
    HPolytope(int d): _d(d) {
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
    MT get_mat() {
        return A;
    }


    // return the vector b
    VT get_vec() {
        return b;
    }


    // change the matrix A
    void set_mat(MT A2) {
        A = A2;
    }


    // change the vector b
    void set_vec(VT b2) {
        b = b2;
    }


    // set a specific coeff of matrix A
    FT get_mat_coeff(int i, int j) {
        return A(i,j);
    }


    // get a spesific coeff of vector b
    FT get_vec_coeff(int i) {
        return b(i);
    }


    // get a specific coeff of matrix A
    void put_mat_coeff(int i, int j, FT value) {
        A(i,j) = value;
    }


    // set a spesific coeff of vector b
    void put_vec_coeff(int i, FT value) {
        b(i) = value;
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
                A(i - 1, j - 1) = -Pin[i][j];
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
        viterator rit, vit;

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
        viterator rit;

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

        viterator lamdait = lamdas.begin(), rit;
        FT lamda = 0, min_plus = FT(maxNT), max_minus = FT(minNT);
        FT sum_nom, sum_denom, c_rand_coord, c_rand_coord_prev;
        int i, j, m = num_of_hyperplanes();

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


    // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    void linear_transformIt(MT T) {
        A = A * T;
    }


    // shift polytope by a point c
    void shift(VT c){
        b = b - A*c;
    }


    // return for each facet the distance from the origin
    std::vector<FT> get_dists(FT radius){
        int i=0;
        std::vector <FT> dists(num_of_hyperplanes(), FT(0));
        typename std::vector<FT>::iterator disit = dists.begin();
        for ( ; disit!=dists.end(); disit++, i++)
            *disit = b(i)/A.row(i).norm();

        return dists;
    }


    // no points given for the rounding, you have to sample from the polytope
    template <class T>
    bool get_points_for_rounding (T &randPoints) {
        return false;
    }
};



//-------------------------------------------------------//
//----------------------V Polytope-----------------------//
//-------------------------------------------------------//


// V-Polytope class
template <typename FT>
class VPolytope{
public:
    typedef Eigen::Matrix<FT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<FT,Eigen::Dynamic,1> VT;

private:
    MT V;  //matrix V. Each row contains a vertex
    VT b;  // vector b that contains first column of ine file
    int _d;  //dimension

public:
    VPolytope() {}

    // return dimension
    int dimension() {
        return _d;
    }


    // this function returns 0. The main sampler requests this function to set the length of lambdas vector
    int num_of_hyperplanes() {
        return 0;
    }


    // compute the number of facets of the cyclic polytope in dimension _d with the same number of vertices
    // this is an upper bound for the number of the facets from McMullen's Upper Bound Theorem
    int upper_bound_of_hyperplanes() {
        int k = num_of_vertices(), d1 = int(std::floor((_d + 1) / 2)), d2 = int(std::floor((_d + 2) / 2));
        int num_of_hyp = 0, nom = 1, denom = 1;
        for (int i = (k - _d + 1); i <= k - d1; i++) {
            nom *= i;
        }
        for (int i = 1; i <= _d - d1; i++) {
            denom *= i;
        }
        num_of_hyp += nom / denom;
        nom = 1;
        denom = 1;

        for (int i = (k - _d + 1); i <= k - d2; i++) {
            nom *= i;
        }
        for (int i = 1; i <= _d - d2; i++) {
            denom *= i;
        }
        num_of_hyp += nom / denom;

        return num_of_hyp;
    }


    // return the number of vertices
    int num_of_vertices() {
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
    FT get_mat_coeff(int i, int j) {
        return V(i,j);
    }


    // get a specific coeff of vector b
    FT get_vec_coeff(int i) {
        return b(i);
    }


    // set a specific coeff of matrix V
    void put_mat_coeff(int i, int j, FT value) {
        V(i,j) = value;
    }


    // set a specific coeff of vector b
    void put_vec_coeff(int i, FT value) {
        b(i) = value;
    }


    // Construct matrix V which contains the vertices row-wise
    void init(std::vector<std::vector<FT> > Pin) {
        _d = Pin[0][1] - 1;
        V.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (int j = 1; j < _d + 1; j++) {
                V(i - 1, j - 1) = Pin[i][j];
            }
        }
    }


    // print polytope in input format
    void print() {
        std::cout << " " << V.rows() << " " << _d << " float" << std::endl;
        for (int i = 0; i < V.rows(); i++) {
            for (int j = 0; j < _d; j++) {
                std::cout << V(i, j) << " ";
            }
            std::cout<<"\n";
        }
    }


    // take d+1 points as input and compute the chebychev ball of the defined simplex
    // done is true when the simplex is full dimensional and false if it is not
    std::pair<Point,FT> get_center_radius_inscribed_simplex(std::vector<Point>::iterator it_beg, std::vector<Point>::iterator it_end, bool &done) {

        Point p0 = *it_beg,p1,c;
        int dim = p0.dimension(),i,j;
        std::vector <FT> temp_p;
        FT radius = 0.0, gi, sum = 0.0;
        MT B(dim,dim);
        MT Bg(dim,dim);
        MT e(1,dim);
        VT row(dim);
        VT g(dim);
        std::pair <Point,FT> result;

        for (j=1; j<dim+1; j++) {
            Point pk = *(it_beg + j);
            e(j - 1) = 1.0;
            for (i = 0; i < dim; i++) {
                B(i, j - 1) = pk[i] - p0[i];
            }
        }
        Bg = B;
        Eigen::FullPivLU <MT> lu_decomp(B);
        auto rank = lu_decomp.rank();
        if(rank==dim) {  // check if the simplex is full dimensional
            done = true;
        }else {
            return result;
        }
        B = B.inverse();
        for (i=0; i<dim; i++) {
            for (j = 0; j < dim; j++) {
                row(j) = B(i, j);
            }
            gi = row.norm();
            radius += gi;
            g(i) = gi;
            if (i < dim - 1) {
                sum += gi;
            }
        }
        e = e * B;
        radius += e.norm();
        radius = 1.0 / radius;
        g = Bg * g;
        g = radius * g;
        for (i=0; i<dim; i++) {
            temp_p.push_back(p0[i] + g(i));
        }
        c = Point(dim, temp_p.begin(), temp_p.end());
        result.first = c;
        result.second = radius;

        return result;
    }


    // pick d+1 random vertices until they define a full dimensional simplex and then compute the chebychev ball of
    // that simplex
    std::pair<Point,FT> chebyshev_center() {

        std::vector<Point> verts(_d+1);
        std::vector<FT> vecp(_d);
        int m = num_of_vertices(), vert_rand, pointer=0,i,j;
        std::vector<int> x_vec(_d);
        bool done=false;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_int_distribution<> uidist(1, m);

        std::pair<Point,FT> res;
        // while d+1 points do not define a full dimensional simplex repeat
        while(!done){
            pointer=0;
            x_vec.assign(_d+1,0);
            while(pointer!=_d+1) {
                vert_rand = uidist(rng);
                // Check if this vertex is selected first time
                if (std::find(x_vec.begin(), x_vec.begin() + pointer, vert_rand) == x_vec.begin() + pointer) {
                    x_vec[pointer] = vert_rand;
                    pointer++;
                }
            }

            for(i=0; i<(_d+1); i++){
                for (j=0; j<_d; j++) {
                    vecp[j] = V(x_vec[i] - 1, j);
                }
                verts[i] = Point(_d,vecp.begin(),vecp.end());
            }
            res=get_center_radius_inscribed_simplex(verts.begin(), verts.end(), done);
        }
        return res;
    }


    // check if point p belongs to the convex hull of V-Polytope P
    int is_in(Point p) {
        if(memLP_Vpoly(V, p)){
            return -1;
        }
        return 0;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with V-polytope
    std::pair<FT,FT> line_intersect(Point r,
                                          Point v) {
        FT min_plus, max_minus;

        max_minus = intersect_line_Vpoly(V, r, v, true);
        min_plus = intersect_line_Vpoly(V, r, v, false);

        return std::pair<FT, FT>(min_plus, max_minus);
    }


    // Compute the intersection of a coordinate ray
    // with a V-polytope
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          int rand_coord,
                                          std::vector<FT> &lamdas) {
        FT min_plus, max_minus;
        std::vector<FT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        max_minus = intersect_line_Vpoly(V, r, v, true);
        min_plus = intersect_line_Vpoly(V, r, v, false);

        return std::pair<FT, FT> (min_plus, max_minus);
    }


    // Compute the intersection of a coordinate ray
    // with a V-polytope
    std::pair<FT,FT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<FT> &lamdas) {
        FT min_plus, max_minus;
        std::vector<FT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        max_minus = intersect_line_Vpoly(V, r, v, true);
        min_plus = intersect_line_Vpoly(V, r, v, false);

        return std::pair<FT, FT> (min_plus, max_minus);
    }


    // shift polytope by a point c
    void shift(VT c) {
        MT V2 = V.transpose().colwise() - c;
        V = V2.transpose();
    }


    // apply linear transformation, of square matrix T, to a V-Polytope
    void linear_transformIt(MT T) {
        MT V2 = T.inverse() * V.transpose();
        V = V2.transpose();
    }


    // consider an upper bound for the number of facets of a V-polytope
    // for each facet consider a lower bound for the distance from the origin
    // useful for CV algorithm to get the first gaussian
    std::vector<FT> get_dists(FT radius) {
        std::vector <FT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }


    // in number_of_vertices<=20*dimension use the vertices for the rounding
    // otherwise you have to sample from the V-polytope
    template <class T>
    bool get_points_for_rounding (T &randPoints) {
        if (num_of_vertices()>20*_d) {
            return false;
        }
        int j;
        std::vector<FT> temp(_d,FT(0));
        typename std::vector<FT>::iterator pointIt;
        for (int i=0; i<num_of_vertices(); i++) {
            pointIt = temp.begin(); j=0;
            for ( ; pointIt!=temp.end(); pointIt++, j++){
                *pointIt = V(i,j);
            }
            Point p(_d, temp.begin(), temp.end());
            randPoints.push_back(p);
        }
        return true;
    }
};


#endif
