//
// Created by Bento Natura on 11/06/2020.
//

#ifndef NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H
#define NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H

//TODO: Separate the main two classes in this file.
#include <vector>
#include "MonomialsClass.h"
#include "NonSymmetricIPM.h"

#include "matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<std::pair<Double, Double> > HyperRectangle;
typedef Vector PolynomialSOS;

class EnvelopeProblemSOS {
public:
    EnvelopeProblemSOS(int num_variables, int max_degree, HyperRectangle &hyperRectangle_);

    //FIXME: Rename as currently the degree of the polynomial remains unchanged.
    static PolynomialSOS prod_sos(PolynomialSOS p1, PolynomialSOS p2) {
        assert(p1.rows() == p2.rows());
        PolynomialSOS p = Vector::Zero(p1.rows());
        for (Eigen::Index i = 0; i < p.rows(); ++i) {
            for (Eigen::Index j = 0; j <= i; j++) {
                p(i) += p1(j) * p2(i - j);
            }
        }
        return p;
    }

    PolynomialSOS generate_zero_polynomial();

    void add_polynomial(PolynomialSOS &polynomial);

    Instance construct_SOS_instance();

    void print_solution(Solution sol);

    void plot_polynomials_and_solution(const Solution &sol);

    Matrix get_transformation_matrix();

private:
    unsigned _n;
    unsigned _d;
    unsigned _L;
    unsigned _U;
    Vector _objectives_vector;
    std::vector<PolynomialSOS> _polynomials_bounds;
    HyperRectangle _hyperRectangle;
    std::vector<PolynomialSOS> _basis_polynomials;
};

#endif //NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H
