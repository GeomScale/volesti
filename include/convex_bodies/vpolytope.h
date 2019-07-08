// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VPOLYTOPE_H
#define VPOLYTOPE_H

#include <limits>

#include <iostream>
#include "solve_lp.h"

//min and max values for the Hit and Run functions

// V-Polytope class
template <class Point, class  RNGType>
class VPolytope{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef RNGType rngtype;

private:
    MT V;  //matrix V. Each row contains a vertex
    VT b;  // vector b that contains first column of ine file
    unsigned int _d;  //dimension
    REAL *conv_comb;
    MT Fmat;

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
        conv_comb = (REAL *) malloc((_d+1) * sizeof(*conv_comb));
        Fmat.resize(V.rows(),_d);
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
        Fmat.resize(_d,_d);
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
    std::pair<Point,NT> get_center_radius_inscribed_simplex(typename std::vector<Point>::iterator it_beg,
                                                            typename std::vector<Point>::iterator it_end, bool &done) {

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
    std::pair<NT,NT> line_intersect(Point r, Point v) {

        return intersect_double_line_Vpoly<NT>(V, r, v);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        return intersect_double_line_Vpoly<NT>(V, r, v);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        return intersect_double_line_Vpoly<NT>(V, r, v);
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        std::pair<NT, int> vppair;
        vppair.first = intersect_line_Vpoly(V, r, v, conv_comb, false, false);
        vppair.second = 1;
        return vppair;
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        std::pair<NT, int> vppair;
        vppair.first = intersect_line_Vpoly(V, r, v, conv_comb, false, false);
        vppair.second = 1;
        return vppair;
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());
        return intersect_double_line_Vpoly<NT>(V, r, v);
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
        return intersect_double_line_Vpoly<NT>(V, r, v);
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

    void compute_reflection(Point &v, Point &p, int facet) {

        int count = 0;
        VT bb = VT::Ones(_d), pp(_d);
        for (int j = 0; j < num_of_vertices(); ++j) {
            if (*(conv_comb + j) > 0.0) {
                //std::cout<<"get vertex "<<*(conv_comb + j)<<std::endl;
                Fmat.row(count) = V.row(j);
                count++;
            } else {
                //std::cout<<"dont get vertex "<<*(conv_comb + j)<<std::endl;
                pp = V.row(j);
            }
        }

        //std::cout<<Fmat<<"\n"<<std::endl;
        VT a = Fmat.colPivHouseholderQr().solve(bb);
        if (a.dot(pp) > 1.0) a = -a;
        a = a/a.norm();
        //std::cout<<"a = "<<a<<"\n"<<std::endl;

        Point s(_d);
        for (int i = 0; i < _d; ++i) {
            s.set_coord(i, a(i));
        }
        s = ((-2.0 * v.dot(s)) * s);
        v = s + v;


    }

};

#endif
