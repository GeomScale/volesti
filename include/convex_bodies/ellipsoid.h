// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis
// Copyright (c) 2021 Vaibhav Thakkar

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
//Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

#include <iostream>
#include <Eigen/Eigen>
#include "volume/math_helpers.hpp"


template <class Point>
class Ellipsoid{
public:
typedef Point PointType;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    // representation is (x - c)' A (x - c) <= 1, center is assumed to be origin for now
    MT A;
    Point c;

    // MT L;   // LL' = A
    unsigned int _dim;

    // eigen vectors and values
    VT _eigen_values;
    VT _eigen_values_inv;
    VT _eigen_values_inv_sqrt;
    MT _Eigen_Vectors;

public:

    Ellipsoid(MT& Ain) : A(Ain) {
        Eigen::SelfAdjointEigenSolver<MT> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error("Eigen solver returned error!");
        }

        _eigen_values = eigensolver.eigenvalues();
        _Eigen_Vectors = eigensolver.eigenvectors();

        _eigen_values_inv = _eigen_values.array().inverse().matrix();
        _eigen_values_inv_sqrt = _eigen_values_inv.array().sqrt().matrix();

        _dim = A.rows();
        c = Point(_dim);
    }


    // Constructor for copula ellipsoid only
    Ellipsoid(std::vector<std::vector<NT> >& Ain) {
        _dim = Ain.size();
        A.resize(_dim, _dim);
        for (unsigned int i = 0; i < Ain.size(); i++) {
            for (unsigned int j = 0; j < Ain.size(); j++) {
                A(i,j) = Ain[i][j];
            }
        }
    }


    VT eigenvalues() const {
        return _eigen_values;
    }


    VT eigenvalues_inv() const {
        return _eigen_values_inv;
    }


    VT eigenvalues_inv_sqrt() const {
        return _eigen_values_inv_sqrt;
    }


    MT eigenvectors() const {
        return _Eigen_Vectors;
    }


    unsigned int dimensions() const {
        return _dim;
    }


    void print() {
        std::cout << "Ellipse is in the form: x' A x <= 1, (center is assumed to be origin always) \n";
        std::cout << "c = \n" << c.print();
        std::cout << "A = \n" << A;
    }



    NT mat_mult(Point const& p) const {
        VT x = p.getCoefficients();
        return x.transpose() * A.template selfadjointView<Eigen::Upper>() * x;
    }


    VT vec_mult(VT const& b) const {
        return A.template selfadjointView<Eigen::Upper>()*b;
    }


    NT log_volume() const {
        NT ball_log_vol = (NT(_dim)/NT(2) * std::log(M_PI)) - log_gamma_function(NT(_dim) / NT(2) + 1);
        NT det_factor = - 0.5 * std::log( A.determinant() );

        return det_factor + ball_log_vol;
    }


    void scale(NT scale_factor) {
        assert (scale_factor > 0);

        NT scale_factor_sq = scale_factor * scale_factor;
        NT inv_scale_factor = (NT(1.0) / scale_factor);
        NT inv_scale_factor_sq = (NT(1.0) / scale_factor_sq);

        // L = mult_factor * L;
        _eigen_values = inv_scale_factor_sq * _eigen_values;
        _eigen_values_inv = scale_factor_sq * _eigen_values_inv;
        _eigen_values_inv_sqrt = scale_factor * _eigen_values_inv_sqrt;

        A = inv_scale_factor_sq * A; // as volume depends on square root of it's determinant
    }


    int is_in(Point const& p) const {
        return mat_mult(p) > 1 ? 0 : -1;
    }


    // compute intersection point of ray starting from r and pointing to v
    std::pair<NT, NT> line_intersect(Point const& r, Point const& v) const {
        // constants of a quadratic equation
        NT a_q = mat_mult(v);
        NT b_q = 2 * r.getCoefficients().dot(vec_mult(v.getCoefficients()));
        NT c_q = mat_mult(r);

        NT D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }


    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    const VT &Ar,
                                    const VT &Av) const
    {
        return line_intersect(r, v);
    }


    std::pair<NT,NT> line_intersect(Point const& r,
                                    Point const& v,
                                    VT& Ar,
                                    VT& Av,
                                    NT &lambda_prev) const
    {
        return line_intersect(r, v);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v) const
    {
        return std::pair<NT,NT>(line_intersect(r, v).first, 0);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              VT& Ar,
                                              VT& Av) const
    {
        return line_positive_intersect(r, v);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              VT& Ar,
                                              VT& Av,
                                              NT const& lambda_prev) const
    {
        return line_positive_intersect(r, v);
    }


    // Compute the intersection of a coordinate ray
    std::pair<NT,NT> line_intersect_coord(Point const& r, const unsigned int rand_coord) const {
        NT a_q = A(rand_coord, rand_coord);
        NT b_q = 2 * r.getCoefficients().dot(A.col(rand_coord));
        NT c_q = mat_mult(r);

        NT D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }


    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          VT& lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }


    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          VT& lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }


    void compute_reflection (Point& v, Point const& p) const
    {
        // normal vector is Ap
        Point s(vec_mult(p.getCoefficients()));
        s *= (1.0 / std::sqrt(s.squared_length()));
        s *= (-2.0 * v.dot(s));
        v += s;
    }
};

#endif
