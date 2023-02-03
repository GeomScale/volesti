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

/// This class represents an ellipsoid parameterized by a point type
/// \tparam Point Point type
template <class Point>
class Ellipsoid{
public:
typedef Point PointType;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

private:
    // representation is x'A x <= 1, i.e center is always assumed to be origin.
    MT A;

    unsigned int _dim;
    MT _L_cov;   // LL' = inv(A) for sampling procedures

    // eigen vectors and values
    VT _eigen_values;
    VT _eigen_values_inv;
    VT _eigen_values_inv_sqrt;
    MT _Eigen_Vectors;

public:

    Ellipsoid() {}

    // TODO(vaithak): Add a flag for telling whether the matrix passed is already inverse
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

        Eigen::LLT<MT> lltOfA(A.inverse()); // compute the Cholesky decomposition of inv(A)
        if (lltOfA.info() != Eigen::Success) {
            throw std::runtime_error("Cholesky decomposition failed!");
        }
        _L_cov = lltOfA.matrixL();
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


    NT radius() const {
        return _eigen_values_inv_sqrt(0);
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


    MT Lcov() const {
        return _L_cov;
    }


    // return L_cov * x
    VT mult_Lcov(VT const& x) const {
        return _L_cov.template triangularView<Eigen::Lower>() * x;
    }


    void print() const {
        std::cout << "Ellipse is in the form: x' A x <= 1, (center is assumed to be origin always) \n";
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
        NT det_factor = std::log( _eigen_values_inv_sqrt.prod() );

        return det_factor + ball_log_vol;
    }


    void scale(NT scale_factor) {
        assert (scale_factor > 0);

        NT scale_factor_sq = scale_factor * scale_factor;
        NT inv_scale_factor = (NT(1.0) / scale_factor);
        NT inv_scale_factor_sq = (NT(1.0) / scale_factor_sq);

        _eigen_values = inv_scale_factor_sq * _eigen_values;
        _eigen_values_inv = scale_factor_sq * _eigen_values_inv;
        _eigen_values_inv_sqrt = scale_factor * _eigen_values_inv_sqrt;
        _L_cov = scale_factor * _L_cov;

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
        NT c_q = mat_mult(r) - 1;

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
        NT res = line_intersect(r, v).first;
        return std::pair<NT,int>(res, 0);
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
        NT c_q = mat_mult(r) - 1;

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
        s *= (1.0 / s.length());
        s *= (-2.0 * v.dot(s));
        v += s;
    }


    template <typename update_parameters>
    void compute_reflection (Point& v, Point const& p, update_parameters &params) const
    {
        // normal vector is Ap
        Point s(vec_mult(p.getCoefficients()));
        params.ball_inner_norm = s.length();

        params.inner_vi_ak = v.dot(s) / params.ball_inner_norm;
        v += (s * (-2.0 * params.inner_vi_ak * (1.0 / params.ball_inner_norm)));
    }
};

#endif