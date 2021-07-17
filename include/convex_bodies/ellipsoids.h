// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

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

    // representation is (x - c)' A (x - c) <= 1
    MT A;
    Point c;

    MT L;   // LL' = A
    unsigned int dim;

public:

    Ellipsoid(MT& Ain) : A(Ain) {
        Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of Ain
        if(lltOfA.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Possibly non semi-positive definitie matrix!");
        }
        L = lltOfA.matrixL();

        dim = A.rows();
        c = Point(dim);
        c.set_to_origin();
    }


    Ellipsoid(MT& Ain, Point& center) : A(Ain), c(center) {
        Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of Ain
        if(lltOfA.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Possibly non semi-positive definitie matrix!");
        }
        L = lltOfA.matrixL();

        dim = A.rows();
    }


    // Constructor for copula ellipsoid
    Ellipsoid(std::vector<std::vector<NT> >& Ain) {
        dim = Ain.size();
        A.resize(dim, dim);
        for (unsigned int i = 0; i < Ain.size(); i++) {
            for (unsigned int j = 0; j < Ain.size(); j++) {
                A(i,j) = Ain[i][j];
            }
        }

        dim = A.rows();
        c = Point(dim);
        c.set_to_origin();
    }


    NT mat_mult(Point const& p) {
        NT res = NT(0);
        for (i=0; i<dim; i++) {
            for (j = i; j<dim; j++) {
                res += 2*A(i, j)*p[i]*p[j];
            }
        }

        return res;
    }

    
    NT log_volume () {
        NT ball_log_vol = (NT(n)/NT(2) * std::log(M_PI)) - log_gamma_function(NT(n) / NT(2) + 1);
        NT det_factor = - std::log( A.determinant() );

        return det_factor + ball_log_vol;
    }


    void scale(NT scale_factor) {
        assert (scale_factor > 0);

        L = (1.0 / scale_factor) * L;
        A = (1.0 / (scale_factor * scale_factor)) * A;
    }


    int is_in(Point const& p){
        NT val = mat_mult(p - c);
        if (val > 1) {
            return 0;
        }

        return -1;
    }


    // compute intersection point of ray starting from r and pointing to v
    std::pair<NT, NT> line_intersect(Point const& r, Point const& v) const
        // constants of a quadratic equation
        NT a_q = mat_mult(v);
        NT b_q = 2 * (r - c).getCoefficients().dot(A * v.getCoefficients());
        NT c_q = mat_mult(r - c);
        
        D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }


    // Compute the intersection of a coordinate ray
    std::pair<NT,NT> line_intersect_coord(const Point &r, const unsigned int rand_coord) const {
        NT a_q = A(rand_coord, rand_coord);
        NT b_q = 2 * (r - c).getCoefficients().dot(A.col(rand_coord));
        NT c_q = mat_mult(r - c);
        
        D = std::pow(b_q, 2) - 4*a_q*c_q;
        return std::pair<NT, NT> ( (-b_q + std::sqrt(D))/(2*a_q) , (-b_q - std::sqrt(D))/(2*a_q) );
    }
};

#endif
