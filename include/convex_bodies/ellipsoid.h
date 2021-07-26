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
#include <boost/math/special_functions/gamma.hpp>


template <typename NT>
NT log_gamma_function(NT x)
{
    if (x <= NT(100)) return std::log(tgamma(x));
    return (std::log(x - NT(1)) + log_gamma_function(x - NT(1)));
}


template <class Point, class MT>
class Ellipsoid{
private:
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    // representation is (x - c)' A (x - c) <= 1, center is assumed to be origin for now
    MT A;
    Point c;

    // MT L;   // LL' = A
    unsigned int _dim;

    // eigen vectors and values
    VT _eigen_values;
    VT _eigen_values_inv;
    VT _eigen_values_inv_sqrt;
    MT _eigen_vecs;

public:

    Ellipsoid(MT& Ain) : A(Ain) {
        // Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of Ain
        // if(lltOfA.info() == Eigen::NumericalIssue) {
        //     throw std::runtime_error("Possibly non semi-positive definitie matrix!");
        // }
        // L = lltOfA.matrixL();

        Eigen::SelfAdjointEigenSolver<MT> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error("Eigen solver returned error!");
        }

        _eigen_values = eigensolver.eigenvalues();
        _eigen_vecs = eigensolver.eigenvectors();

        _eigen_values_inv = _eigen_values.array().inverse().matrix();
        _eigen_values_inv_sqrt = _eigen_values_inv.array().sqrt().matrix();

        _dim = A.rows();
        c = Point(_dim);
    }


    // Constructor for copula ellipsoid
    Ellipsoid(std::vector<std::vector<NT> >& Ain) {
        _dim = Ain.size();
        A.resize(_dim, _dim);
        for (unsigned int i = 0; i < Ain.size(); i++) {
            for (unsigned int j = 0; j < Ain.size(); j++) {
                A(i,j) = Ain[i][j];
            }
        }

        Eigen::SelfAdjointEigenSolver<MT> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) {
            throw std::runtime_error("Eigen solver returned error!");
        }

        _eigen_values = eigensolver.eigenvalues();
        _eigen_vecs = eigensolver.eigenvectors();

        _eigen_values_inv = _eigen_values.array().inverse().matrix();
        _eigen_values_inv_sqrt = _eigen_values_inv.array().sqrt().matrix();

        _dim = A.rows();
        c = Point(_dim);
    }


    VT eigenvals() const {
        return _eigen_values;
    }


    VT eigenvals_inv() const {
        return _eigen_values_inv;
    }


    VT eigenvals_inv_sqrt() const {
        return _eigen_values_inv_sqrt;
    }


    MT eigenvecs() const {
        return _eigen_vecs;
    }


    unsigned int dimensions() const {
        return _dim;
    }


    void print() {
        std::cout << "Ellipse is in the form: x' A x <= 1, (center is assumed to be origin always) \n";
        std::cout << "c = \n" << c.print();
        std::cout << "A = \n" << A;
    }



    NT mat_mult(Point const& p) {
        VT x = p.getCoefficients();

        if (_dim < 15) {
            NT sum = 0;
            for (Eigen::Index i = 0; i < _dim; ++i) {
                const auto x_i = x[i];
                sum += A.coeff(i, i) * x_i * x_i;
                for (Eigen::Index j = 0; j < i; ++j) {
                    sum += 2 * A.coeff(j, i) * x_i * x[j];
                }
            }
            return sum;
        }

        return x.transpose() * A.template selfadjointView<Eigen::Upper>() * x;
    }


    VT vec_mult(VT& b) {
        return A.template selfadjointView<Eigen::Upper>()*b;
    }


    NT log_volume () {
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


    int is_in(Point const& p){
        NT val = mat_mult(p);
        if (val > 1) {
            return 0;
        }

        return -1;
    }


    // compute intersection point of ray starting from r and pointing to v
    std::pair<NT, NT> line_intersect(Point& r, Point& v) const {
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
                                    const VT &Ar,
                                    const VT &Av,
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
                                              const VT &Ar,
                                              const VT &Av) const
    {
        return line_positive_intersect(r, v);
    }


    std::pair<NT,int> line_positive_intersect(Point const& r,
                                              Point const& v,
                                              const VT &Ar,
                                              const VT &Av,
                                              NT &lambda_prev) const
    {
        return line_positive_intersect(r, v);
    }


    // Compute the intersection of a coordinate ray
    std::pair<NT,NT> line_intersect_coord(const Point &r, const unsigned int rand_coord) const {
        NT a_q = A(rand_coord, rand_coord);
        NT b_q = 2 * r.getCoefficients().dot(A.col(rand_coord));
        NT c_q = mat_mult(r);

        NT D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }


    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          unsigned int const& rand_coord,
                                          const VT &lamdas) const
    {
        return line_intersect_coord(r, rand_coord);
    }


    std::pair<NT,NT> line_intersect_coord(Point const& r,
                                          Point const& r_prev,
                                          unsigned int const& rand_coord,
                                          unsigned int const& rand_coord_prev,
                                          const VT &lamdas) const
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
