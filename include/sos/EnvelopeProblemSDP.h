// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSDP_H
#define NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSDP_H

//TODO: Separate the main two classes in this file.
#include <vector>
#include "MonomialsClass.h"
#include "NonSymmetricIPM.h"

#include "matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<std::pair<IPMDouble, IPMDouble> > HyperRectangle;
typedef Matrix PolynomialSDP;


class EnvelopeProblemSDP {
public:

    EnvelopeProblemSDP(int num_variables, int max_degree, HyperRectangle &hyperRectangle_);

    Matrix get_objective_matrix() const;

    void add_polynomial(PolynomialSDP &polynomial);

    PolynomialSDP generate_zero_polynomial();

    void construct_polynomial_matrix();

    //Requires that polynomial matrix was already constructed
    void construct_objective_matrix();

    IPMDouble calculate_objective(Monomial m);

    IPMDouble calculate_objective(Monomial m, unsigned var);

    //FIXME: Remove trivial rows. Also, the variables Y might not be necessary. The whole barrier can be applied to X itself;
    Instance construct_SDP_instance();

    void print_solution(const Solution &sol);

    Monomial get_matrix_entry(unsigned const row, unsigned const col);

    void print_polynomial(Matrix M) const;

    //Assume that monomial is ordered as 1, x, x*x, ...
    IPMDouble univariate_monomial_evaluation(Monomial const m, IPMDouble const x);

    IPMDouble univariate_polynomial_evaluation(PolynomialSDP const poly, IPMDouble x);

    void plot_polynomials_and_solution(const Solution &sol);

    MonomialsClass monomialObject;

private:
    unsigned _n;
    unsigned _d;
    std::vector<std::vector<Monomial> > _polynomial_product_matrix;
    Matrix _objectives_matrix;
    std::vector<PolynomialSDP> _polynomials_bounds;
    HyperRectangle _hyperRectangle;
};

#endif //NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSDP_H
