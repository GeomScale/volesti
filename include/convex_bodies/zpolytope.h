// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ZPOLYTOPE_H
#define ZPOLYTOPE_H

#include <limits>

#include <iostream>
#include <Eigen/Eigen>
#include "lp_oracles/vpolyoracles.h"
#include "lp_oracles/zpolyoracles.h"

template <typename Point>
class Zonotope {
public:
    typedef Point PointType;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

private:
    MT V;  //matrix V. Each row contains a vertex
    VT b;  // vector b that contains first column of ine file
    unsigned int _d;  //dimension
    std::pair<Point,NT> _inner_ball;
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();

    REAL *conv_comb, *row, *row_mem;
    int *colno, *colno_mem;
    MT sigma;
    MT Q0;
    MT T;

public:

    Zonotope() {}

    // return the dimension
    unsigned int dimension() const
    {
        return _d;
    }

    // this function returns 0. The main sampler requests this function to set the length of lambdas vector
    int num_of_hyperplanes() const
    {
        return 0;
    }


    // return the number of parallelopipeds. Used in get_dists fnction.
    unsigned int upper_bound_of_hyperplanes() const
    {
        return 2*_d;
    }

    void compute_eigenvectors(MT const& G)
    {

        int k = G.cols();
        MT ps = G;
        sigma.resize(k,k);
        sigma = ps.transpose()*ps;
        sigma = (sigma + sigma.transpose()) / 2;
        Eigen::SelfAdjointEigenSolver<MT> es(sigma);

        MT D = es.eigenvalues().asDiagonal();
        MT Q2 = es.eigenvectors();

        Q0.resize(k,k-_d);
        int count=0;
        for (int i = 0; i < k; ++i)
        {
            if (es.eigenvalues()[i]<0.0000001)
            {
                for (int j = 0; j < k; ++j)
                {
                    Q0(j, count) = Q2(j, i);
                }
                count++;
            }
        }
        Eigen::JacobiSVD<MT> svd(Q0, Eigen::ComputeFullU | Eigen::ComputeFullV);
        MT T2 = svd.matrixU().transpose();
        T.resize(_d,k);
        for (int i = k-_d; i < k; ++i)
        {
            for (int j = 0; j < k; ++j)
            {
                T(i-k+_d,j) = T2(i,j);
            }
        }

        for (int i1 = 0; i1 < k; ++i1)
        {
            sigma(i1,i1) = sigma(i1,i1) + 0.00000001;
        }
    }

    MT get_T() const
    {
        return T;
    }

    MT get_Q0() const
    {
        return Q0;
    }

    MT get_sigma() const
    {
        return sigma;
    }

    // return the number of vertices
    int num_of_vertices() const
    {
        return V.rows();
    }

    std::pair<Point,NT> InnerBall() const
    {
        return _inner_ball;
    }

    // return the number of generators
    int num_of_generators() const
    {
        return V.rows();
    }


    // return the matrix V
    MT get_mat() const
    {
        return V;
    }


    // return the vector b
    VT get_vec() const
    {
        return b;
    }


    // change the matrix V
    void set_mat(MT const& V2)
    {
        V = V2;
    }


    // change the vector b
    void set_vec(VT const& b2)
    {
        b = b2;
    }

    Point get_mean_of_vertices() const
    {
        return Point(_d);
    }


    NT get_max_vert_norm() const
    {
        return 0.0;
    }

    // define zonotope using Eigen matrix V. Vector b is neded in order
    // the code to compatible with Hpolytope class
    void init(unsigned int const dim, MT const& _V, VT const& _b)
    {
        _d = dim;
        V = _V;
        b = _b;
        conv_comb = (REAL *) malloc((V.rows()+1) * sizeof(*conv_comb));
        colno = (int *) malloc((V.rows()+1) * sizeof(*colno));
        row = (REAL *) malloc((V.rows()+1) * sizeof(*row));
        colno_mem = (int *) malloc((V.rows()) * sizeof(*colno_mem));
        row_mem = (REAL *) malloc((V.rows()) * sizeof(*row_mem));
        compute_eigenvectors(V.transpose());
    }


    // Construct matrix V which contains the vertices row-wise
    void init(std::vector<std::vector<NT> > const& Pin)
    {
        _d = Pin[0][1] - 1;
        V.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++)
        {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++)
            {
                V(i - 1, j - 1) = Pin[i][j];
            }
        }
        conv_comb = (REAL *) malloc(Pin.size() * sizeof(*conv_comb));
        colno = (int *) malloc((V.rows()+1) * sizeof(*colno));
        row = (REAL *) malloc((V.rows()+1) * sizeof(*row));
        colno_mem = (int *) malloc((V.rows()) * sizeof(*colno_mem));
        row_mem = (REAL *) malloc((V.rows()) * sizeof(*row_mem));
        compute_eigenvectors(V.transpose());
    }


    // print polytope in input format
    void print()
    {
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
    int is_in(Point const& p) const
    {
        if(memLP_Zonotope(V, p, row_mem, colno_mem))
        {
            return -1;
        }
        return 0;
    }


    // Compute an inner ball of the zonotope
    std::pair<Point,NT> ComputeInnerBall()
    {
        std::vector<NT> temp(_d,0);
        NT radius =  maxNT, min_plus;
        Point center(_d);

        for (unsigned int i = 0; i < _d; ++i) {
            temp.assign(_d,0);
            temp[i] = 1.0;
            Point v(_d,temp.begin(), temp.end());
            min_plus = intersect_line_Vpoly<NT>(V, center, v, conv_comb,
                                                row, colno, false, true);
            if (min_plus < radius) radius = min_plus;
        }

        radius = radius / std::sqrt(NT(_d));
        _inner_ball = std::pair<Point, NT> (center, radius);
        return _inner_ball;
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point const& r, Point const& v) const
    {
        return intersect_line_zono(V, r, v, conv_comb, colno);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT const& Ar,
                                    VT const& Av) const
    {
        return intersect_line_zono(V, r, v, conv_comb, colno);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the Zonotope
    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT const& Ar,
                                    VT const& Av,
                                    NT const& lambda_prev) const
    {
        return intersect_line_zono(V, r, v, conv_comb, colno);
    }

    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT const& Ar,
                                               VT const& Av) const
    {
        return std::pair<NT, int> (intersect_line_Vpoly(V, r, v, conv_comb,
                                                        row, colno,
                                                        false, true), 1);
    }


    std::pair<NT, int> line_positive_intersect(Point const& r,
                                               Point const& v,
                                               VT const& Ar,
                                               VT const& Av,
                                               NT const& lambda_prev) const
    {
        return line_positive_intersect(r, v, Ar, Av);
    }


    // Compute the intersection of a coordinate ray
    // with the Zonotope
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          const unsigned int rand_coord,
                                          VT const& lamdas) const
    {

        std::vector<NT> temp(_d,0);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());

        return intersect_line_zono(V, r, v, conv_comb, colno);

    }


    // Compute the intersection of a coordinate ray
    // with the Zonotope
    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          const unsigned int rand_coord,
                                          const unsigned int rand_coord_prev,
                                          VT const& lamdas) const
    {
        return line_intersect_coord(r, rand_coord, lamdas);
    }


    // shift polytope by a point c
    // vector c has to be always the zero vector
    void shift(VT const& c)
    {
        return;
    }


    // get number of parallelopipeds
    // for each parallelopiped consider a lower bound for the distance from the origin
    // useful for CG algorithm to get the first gaussian
    std::vector<NT> get_dists(NT const& radius) const
    {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    // apply linear transformation, of square matrix T, to thr Zonotope
    void linear_transformIt(MT const& T)
    {
        MT V2 = T.inverse() * V.transpose();
        V = V2.transpose();
    }

    // return false to the rounding function
    // no points are given so they have o be sampled
    template <typename T>
    bool get_points_for_rounding (T const& randPoints)
    {
        return false;
    }

    void normalize() {}

    void compute_reflection(Point &v, Point const& p, int const& facet) const
    {
        int count = 0;
        MT Fmat(_d-1,_d);
        const NT e = 0.0000000001;
        for (int j = 0; j < num_of_generators(); ++j)
        {
            if (((1.0 - *(conv_comb + j) ) > e || (1.0 - *(conv_comb + j) )
                                           > e*std::abs(*(conv_comb + j))) &&
                ((1.0 + *(conv_comb + j) ) > e || (1.0 + *(conv_comb + j) )
                                           > e*std::abs(*(conv_comb + j))))
            {
                Fmat.row(count) = V.row(j);
                count++;
            }
        }

        VT a = Fmat.fullPivLu().kernel();

        if(p.getCoefficients().dot(a) < 0.0) a *= -1.0;

        a = a/a.norm();

        // compute reflection
        a *= (-2.0 * v.dot(a));
        v += a;
    }

    void free_them_all() {
        free(row);
        free(colno);
        free(conv_comb);
        free(row_mem);
        free(colno_mem);
    }

};

#endif
