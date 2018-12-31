// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "vars.h"
#include "solve_lp.h"

#ifdef USE_FAISS
#include "faiss/IndexFlat.h"
#endif

//min and max values for the Hit and Run functions


// H-polytope class
template <class Point>
class HPolytope{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    typedef RowMatrixXd MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    MT A; //matrix A
    VT b; // vector b, s.t.: Ax<=b
    unsigned int            _d; //dimension
    NT maxNT = 1.79769e+308;
    NT minNT = -1.79769e+308;

    NT** representatives;
    double maxDistToBoundary;
#ifdef USE_FAISS
    faiss::IndexFlatL2* index; // call constructor
#endif

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

    HPolytope(Eigen::Map<MT> A, Eigen::Map<VT> b) {
        this->A = A;
        this->b = b;
        _d = A.cols();
        representatives = NULL;
    }

    ~HPolytope() {
        if (representatives!=NULL) {
            for (int i=0; i<num_of_hyperplanes()+1; i++) {
                delete representatives[i];
            }
            delete []representatives;
#ifdef USE_FAISS
            delete index;
#endif
        }

    }

    NT* project(NT* point, int facet_idx) {
        /**
		 * for a hyperplane H:=ax=b and a point p:
         * proj_H(p)=p - a * ((p.dot_product(a) - b)/ a.dot_product(a))
         */
        NT dot_product = 0;
        NT normalizer = 0;
        for (uint i=0; i<dimension(); i++) {
            dot_product += A(facet_idx, i)*point[i];
            normalizer += A(facet_idx, i)*A(facet_idx, i);
        }

        double dir_coeff = (dot_product - b[facet_idx]) / normalizer;

        NT* projected_point = new NT[dimension()];
        for (uint i=0; i<dimension(); i++) {
            projected_point[i] = point[i] - dir_coeff*A(facet_idx, i);
        }
        return projected_point;
    }

    NT* get_reflexive_point(NT* internalPoint, int facet_idx) {
		/**
		 * Get the reflexive point of a point p inside the polytope about the facet_idx-th facet
		 * ie. if the facet F is defined by the supporting hyperplane H := ax=b, then
		 * refl_F(p) = p+2*(proj_H(p)-p) <=> 2*proj_H(p)-p
		 */
        NT* projection = project(internalPoint, facet_idx);
        for (uint i=0; i<dimension(); i++) {
            projection[i] = 2*projection[i] - internalPoint[i];
        }

        return projection;
    }

    double squared_distance(NT* a, NT* b) {
        double dist = 0.0;
        for (uint i=0; i<dimension(); i++) {
            double tmp = a[i]-b[i];
            dist += tmp*tmp;
        }
        return dist;
    }

    void create_point_representation(NT* internalPoint) {
        representatives = new NT*[num_of_hyperplanes()+1];
        maxDistToBoundary=-1;

        representatives[0] = this->get_reflexive_point(internalPoint, 0);
        maxDistToBoundary = squared_distance(internalPoint, representatives[0]);

        for (int i=1; i<num_of_hyperplanes(); i++) {
            representatives[i] = this->get_reflexive_point(internalPoint, i);
            double tmpDist = squared_distance(internalPoint, representatives[i]);
            if (tmpDist>maxDistToBoundary) {
                maxDistToBoundary = tmpDist;
            }
        }

        maxDistToBoundary = std::sqrt(maxDistToBoundary);

        representatives[num_of_hyperplanes()] = new NT[dimension()];
        for (uint i=0; i<dimension(); i++) {
            representatives[num_of_hyperplanes()][i] = internalPoint[i];
        }
#ifdef USE_FAISS
        index = new faiss::IndexFlatL2(dimension());
        for (int i=0;i<num_of_hyperplanes()+1; i++)
            index->add(1, representatives[i]);
#endif
    }

#ifdef USE_FAISS
    long get_nearest_facet(NT* point) {
        long *I = new long[1];
        float *D = new float[1];
        index->search(1, point, 1, D, I);

        long nnIndex = I[0];

        delete []I;
        delete []D;
        return nnIndex;
    }
    bool contains_point(NT* point) {
        long nnIndex = get_nearest_facet(point); 
        return nnIndex==num_of_hyperplanes();
    }
#endif

    Point intersect_ray_hyperplane(Point& source, Point& direction, int facet_idx) {
        //l = (b -a*s)/(r*a);
        double source_dot = 0.0;
        double dir_dot = 0.0;

        for (uint i=0; i<dimension(); i++) {
            source_dot += A(facet_idx,i)*source[i];
            dir_dot += A(facet_idx,i)*direction[i];
        }

        double lambda = (b[facet_idx]-source_dot)/dir_dot;

        Point intersection(dimension());
        for (uint i=0; i<dimension(); i++) {
            intersection.set_coord(i, source[i]+lambda*direction[i]);
        }
        return intersection;
    }

#ifdef asd
    Point compute_boundary_intersection(Point& source, Point& direction, float epsilon, int maxSteps) {
        direction.normalize();
        Point direction_epsilon(dimension());
        for (uint j=0; j<dimension(); j++) {
            direction_epsilon.set_coord(j, direction[j]*(-epsilon));
        }
        Point x0(dimension());
        for (uint i=0; i<dimension(); i++) {
            x0.set_coord(i, source[i]+direction[i]*2*maxDistToBoundary);
        }

        for (int currentIt=0; currentIt<maxSteps; currentIt++) {
			/** Find nearest facet */
            long nn_index = get_nearest_facet(x0.data());

            /** if point is not inside */
            if (nn_index!=num_of_hyperplanes()) {
                Point x1 = intersect_ray_hyperplane(source, direction, nn_index);

                double x1_ray_norm = 0.0;
                double x0_ray_norm = 0.0;
                for (uint d_idx=0; d_idx<dimension(); ++d_idx) {
                    double x1_tmp = x1[d_idx]-source[d_idx];
                    x1_ray_norm += x1_tmp*x1_tmp;

                    double x0_tmp = x0[d_idx]-source[d_idx];
                    x0_ray_norm += x0_tmp*x0_tmp;
                }

                if (x1_ray_norm>=x0_ray_norm) {
                    for (uint j=0; j<dimension(); j++) {
                        x0.set_coord(j, x0[j]+direction_epsilon[j]);
                    }
                }
                else {
                    for (uint j=0; j<dimension(); j++) {
                        x0.set_coord(j, x1[j]);
                    }
                }
            }

            else{
                break;
            }
        }
        return x0;
    }
#endif
    //Check if Point p is in H-polytope P:= Ax<=b
    int is_in(Point p) {
        NT sum;
        int m = A.rows();
        for (int i = 0; i < m; i++) {
            sum = b(i);
            for (unsigned int j = 0; j < _d; j++) {
                sum -= A(i, j) * p[j];
            }
            if (sum < NT(0)) { //Check if corresponding hyperplane is violated
                return 0;
            }
        }
        return -1;
    }

    // return dimension
    unsigned int dimension() {
        return _d;
    }


    // return the number of facets
    int num_of_hyperplanes() {
        return A.rows();
    }


    // return the matrix A
    MT& get_mat() {
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
    NT get_mat_coeff(unsigned int i, unsigned int j) {
        return A(i,j);
    }


    // get a spesific coeff of vector b
    NT get_vec_coeff(unsigned int i) {
        return b(i);
    }


    // get a specific coeff of matrix A
    void put_mat_coeff(unsigned int i, unsigned int j, NT value) {
        A(i,j) = value;
    }


    // set a spesific coeff of vector b
    void put_vec_coeff(unsigned int i, NT value) {
        b(i) = value;
    }


    void init(unsigned int dim, MT _A, VT _b) {
        _d = dim;
        A = _A;
        b = _b;
    }


    // default initialize: cube(d)
    void init(unsigned int d) {
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


    //define matrix A and vector b, s.t. Ax<=b and the dimension
    void init(std::vector<std::vector<NT> > Pin) {
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
                std::cout << A(i, j) << " ";
                #endif
            }
            #ifdef VOLESTI_DEBUG
            std::cout << "<= " << b(i) << std::endl;
            #endif
        }

#ifdef VOLESTI_DEBUG
        for (int i=0; i<num_of_hyperplanes()+1; i++) {
            std::cout<<"#"<<i<<": "<<representatives[i][0];
            for (uint j=1; j<dimension(); j++) {
                std::cout<<", "<<representatives[i][j];
            }
            std::cout<<std::endl;
        }
#endif
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

    

    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point,NT> ComputeInnerBall() {

        std::pair <Point,NT> res;
        res = ComputeChebychevBall<NT, Point>(A, b);  //lpSolve lib for the linear program
        return res;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<Point, Point> line_intersect(Point& r,
                                          Point& v) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes();
        viterator rit, vit;

        for (int i = 0; i < m; i++) {
            sum_nom = b(i);
            sum_denom = NT(0);
            for (uint j=0; j<dimension(); j++){
                sum_nom -= A(i, j) * r[j];
                sum_denom += A(i, j) * v[j];
            }
            if (sum_denom == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                //std::cout << sum_nom << "\t"<<sum_denom<<std::endl;
                lamda = sum_nom / sum_denom;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
        }
        //r+(min_plus*v),r+(max_minus*v)
        Point a(dimension());
        Point b(dimension());

        //std::cout<<"min plus: " <<min_plus<<", max minus: "<<max_minus<<std::endl;
        for (uint i=0; i<dimension(); i++) {
            a.set_coord(i, r[i] + min_plus*v[i]);
            b.set_coord(i, r[i] + max_minus*v[i]);
        }
        return std::pair<Point, Point>(a, b);
    }


    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom;
        unsigned int j;
        int m = num_of_hyperplanes();
        viterator rit;

        for (int i = 0; i < m; i++) {
            sum_nom = b(i);
            sum_denom = A(i, rand_coord);
            rit = r.iter_begin();
            j = 0;
            for (; rit != r.iter_end(); rit++, j++) {
                sum_nom -= A(i, j) * (*rit);
            }
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
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        viterator lamdait = lamdas.begin(), rit;
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        NT sum_nom, sum_denom, c_rand_coord, c_rand_coord_prev;
        int m = num_of_hyperplanes();

        for (int i = 0; i < m; i++) {
            sum_denom = b(i);
            c_rand_coord = A(i, rand_coord);
            c_rand_coord_prev = A(i, rand_coord_prev);

            *lamdait = *lamdait
                       + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
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
    void linear_transformIt(MT T) {
        A = A * T;
    }


    // shift polytope by a point c
    void shift(VT c){
        b = b - A*c;
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(NT radius){
        unsigned int i=0;
        std::vector <NT> dists(num_of_hyperplanes(), NT(0));
        typename std::vector<NT>::iterator disit = dists.begin();
        for ( ; disit!=dists.end(); disit++, i++)
            *disit = b(i) / A.row(i).norm();

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
template <class Point, class  RNGType>
class VPolytope{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    MT V;  //matrix V. Each row contains a vertex
    VT b;  // vector b that contains first column of ine file
    unsigned int _d;  //dimension
    NT maxNT = 1.79769e+308;
    NT minNT = -1.79769e+308;

public:
    VPolytope() {}

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
        unsigned int k = num_of_vertices(), d1 = std::floor((_d + 1) / 2), d2 = std::floor((_d + 2) / 2), res;
        long double num_of_hyp = 0.0, nom = 1.0, denom = 1.0;
        for (unsigned int i = (k - _d + 1); i <= k - d1; i++) {
            nom *= i;
        }
        for (unsigned int i = 1; i <= _d - d1; i++) {
            denom *= i;
        }
        num_of_hyp += nom / denom;
        nom = 1.0;
        denom = 1.0;

        for (unsigned int i = (k - _d + 1); i <= k - d2; i++) {
            nom *= i;
        }
        for (unsigned int i = 1; i <= _d - d2; i++) {
            denom *= i;
        }
        num_of_hyp += nom / denom;

        res = num_of_hyp;
        return res;
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


    void init(unsigned int dim, MT _V, VT _b) {
        _d = dim;
        V = _V;
        b = _b;
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


    // take d+1 points as input and compute the chebychev ball of the defined simplex
    // done is true when the simplex is full dimensional and false if it is not
    std::pair<Point,NT> get_center_radius_inscribed_simplex(typename std::vector<Point>::iterator it_beg, typename std::vector<Point>::iterator it_end, bool &done) {

        Point p0 = *it_beg,p1,c;
        unsigned int dim = p0.dimension();
        unsigned int i,j;
        std::vector <NT> temp_p;
        NT radius = 0.0, gi, sum = 0.0;
        MT B(dim,dim);
        MT Bg(dim,dim);
        MT e(1,dim);
        VT row(dim);
        VT g(dim);
        std::pair <Point,NT> result;

        for (j=1; j<dim+1; j++) {
            Point pk = *(it_beg + j);
            e(j - 1) = 1.0;
            for (i = 0; i < dim; i++) {
                B(i, j - 1) = pk[i] - p0[i];
            }
        }
        Bg = B;
        Eigen::FullPivLU <MT> lu_decomp(B);
        int rank = lu_decomp.rank();
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


    // pick d+1 random vertices until they define a full dimensional simplex and then
    // compute the chebychev ball of that simplex
    std::pair<Point,NT> ComputeInnerBall() {

        std::vector<Point> verts(_d+1);
        std::vector<NT> vecp(_d);
        unsigned int vert_rand, pointer=0;
        unsigned int i,j;
        int m = num_of_vertices();
        std::vector<int> x_vec(_d);
        bool done=false;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::random::uniform_int_distribution<> uidist(1, m);

        std::pair<Point,NT> res;
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
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r,
                                          Point v) {
        NT min_plus, max_minus;

        max_minus = intersect_line_Vpoly<NT>(V, r, v, true, false);
        min_plus = intersect_line_Vpoly<NT>(V, r, v, false, false);

        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {
        NT min_plus, max_minus;
        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        max_minus = intersect_line_Vpoly<NT>(V, r, v, true, false);
        min_plus = intersect_line_Vpoly<NT>(V, r, v, false, false);

        return std::pair<NT, NT> (min_plus, max_minus);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {
        NT min_plus, max_minus;
        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        max_minus = intersect_line_Vpoly<NT>(V, r, v, true, false);
        min_plus = intersect_line_Vpoly<NT>(V, r, v, false, false);

        return std::pair<NT, NT> (min_plus, max_minus);
    }


    // shift polytope by a point c
    void shift(VT c) {
        MT V2 = V.transpose().colwise() - c;
        V = V2.transpose();
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(MT T) {
        MT V2 = T.inverse() * V.transpose();
        V = V2.transpose();
    }


    // consider an upper bound for the number of facets of a V-polytope
    // for each facet consider a lower bound for the distance from the origin
    // useful for CV algorithm to get the first gaussian
    std::vector<NT> get_dists(NT radius) {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }


    // in number_of_vertices<=20*dimension use the vertices for the rounding
    // otherwise you have to sample from the V-polytope
    template <class PointList>
    bool get_points_for_rounding (PointList &randPoints) {
        if (num_of_vertices()>20*_d) {
            return false;
        }
        unsigned int j;
        std::vector<NT> temp(_d,NT(0));
        typename std::vector<NT>::iterator pointIt;
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


//--------------------------//
//-----Zonotope class-------//
//--------------------------//

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
    NT maxNT = 1.79769e+308;
    NT minNT = -1.79769e+308;

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

        res = num_of_hyp;
        return res;
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
        initial_shifting(); // shift zonotope to the origin
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
        initial_shifting(); // shift zonotope to the origin
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
            min_plus = intersect_line_Vpoly<NT>(V, center, v, false, true);
            if (min_plus < radius) radius = min_plus;
        }

        radius = radius / std::sqrt(NT(_d));
        return std::pair<Point, NT> (center, radius);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point r,
                                    Point v) {
        NT min_plus, max_minus;

        max_minus = intersect_line_Vpoly<NT>(V, r, v, true, true);
        min_plus = intersect_line_Vpoly<NT>(V, r, v, false, true);

        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // Compute the intersection of a coordinate ray
    // with the Zonotope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {
        NT min_plus, max_minus;
        std::vector<NT> temp(_d,0);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        max_minus = intersect_line_Vpoly<NT>(V, r, v, true, true);
        min_plus = intersect_line_Vpoly<NT>(V, r, v, false, true);

        return std::pair<NT, NT> (min_plus, max_minus);
    }


    // Compute the intersection of a coordinate ray
    // with the Zonotope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {
        NT min_plus, max_minus;
        std::vector<NT> temp(_d,0);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        max_minus = intersect_line_Vpoly<NT>(V, r, v, true, true);
        min_plus = intersect_line_Vpoly<NT>(V, r, v, false, true);

        return std::pair<NT, NT> (min_plus, max_minus);
    }


    // shift zonotope to the origin
    void initial_shifting() {
        std::vector<NT> vec(_d,0.0);
        typename std::vector<NT>::iterator vecit;
        Point xc(_d);

        int m = V.rows(), j;
        Point temp;

        for (int i = 0; i < m; ++i) {
            vecit = vec.begin();
            j = 0;
            for ( ; vecit!=vec.end(); ++vecit, ++j) {
                *vecit = V(i,j);
            }
            temp = Point(_d, vec.begin(), vec.end());
            xc = xc + temp;
        }
        xc = xc * (1.0 / NT(m));

        VT c(_d);
        vecit = xc.iter_begin();
        j = 0;
        for ( ; vecit!=xc.iter_end(); ++vecit, ++j) {
            c(j) = *vecit;
        }

        MT V2 = V.transpose().colwise() - c;
        V = V2.transpose();
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

};


#endif
