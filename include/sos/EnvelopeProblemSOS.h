// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H
#define NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H

//TODO: Separate the main two classes in this file.
#include <vector>
#include "MonomialsClass.h"
#include "NonSymmetricIPM.h"

#include "matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<std::pair<IPMDouble, IPMDouble> > HyperRectangle;
typedef Vector PolynomialSOS;

class EnvelopeProblemSOS {
public:
    EnvelopeProblemSOS(unsigned num_variables, unsigned max_degree, HyperRectangle &hyperRectangle_);

    //FIXME: Rename as currently the degree of the polynomial remains unchanged.
    static InterpolantVector polynomial_product(InterpolantVector p1, InterpolantVector p2) {
        assert(p1.rows() == p2.rows());
        InterpolantVector p = InterpolantVector::Zero(p1.rows());
        for (Eigen::Index i = 0; i < p.rows(); ++i) {
            for (Eigen::Index j = 0; j <= i; j++) {
                p(i) += p1(j) * p2(i - j);
            }
        }
        return p;
    }

    static InterpolantVector polynomial_product(std::vector<InterpolantVector> const poly_vec) {
        auto len = poly_vec.size();
        assert(not poly_vec.empty());
        if (len == 1) {
            return poly_vec[0];
        }
        if (len == 2) {
            return polynomial_product(poly_vec[0], poly_vec[1]);
        }
        auto mid = len / 2;
        auto first_it = poly_vec.begin();
        auto mid_it = poly_vec.begin() + mid;
        auto last_it = poly_vec.end();
        std::vector<InterpolantVector> vec1(first_it, mid_it);
        std::vector<InterpolantVector> vec2(mid_it, last_it);
        return polynomial_product(polynomial_product(vec1), polynomial_product(vec2));
    }

    InterpolantVector generate_zero_polynomial();

    void add_polynomial(InterpolantVector &polynomial);

    Instance construct_SOS_instance();

    void print_solution(Solution sol);

    void plot_polynomials_and_solution(const Solution &sol);

    void calculate_basis_polynomials();

    InterpolantMatrix get_transformation_matrix();

private:
    unsigned _n;
    unsigned _d;
    unsigned _L;
    unsigned _U;
    InterpolantVector _objectives_vector;
    std::vector<InterpolantVector> _polynomials_bounds;
    HyperRectangle _hyperRectangle;
    std::vector<InterpolantVector> _basis_polynomials;
    std::shared_ptr<spdlog::logger> _logger;

    //This bool sets whether the transformation matrix from the
    //standard monomial basis to the Lagrange basis through the
    //Chebyshev points of second kind is calculated.
    //Set to true if you want to save runtime but have
    //arbitrary polynomials plotted.
    bool _input_in_interpolant_basis = true;

    bool _use_weighted_polynomials = true;

    void get_clenshaw_curtis_integrals();

    void initialize_loggers();
};

#endif //NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H
