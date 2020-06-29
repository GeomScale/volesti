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

typedef std::vector<std::pair<IPMDouble, IPMDouble> > HyperRectangle;
typedef Vector PolynomialSOS;

class EnvelopeProblemSOS {
public:
    EnvelopeProblemSOS(int num_variables, int max_degree, HyperRectangle &hyperRectangle_);

    //FIXME: Rename as currently the degree of the polynomial remains unchanged.
    static InterpolantVector prod_sos(InterpolantVector p1, InterpolantVector p2) {
        assert(p1.rows() == p2.rows());
        InterpolantVector p = InterpolantVector::Zero(p1.rows());
        for (Eigen::Index i = 0; i < p.rows(); ++i) {
            for (Eigen::Index j = 0; j <= i; j++) {
                p(i) += p1(j) * p2(i - j);
            }
        }
        return p;
    }

    static InterpolantVector  prod_sos(std::vector<InterpolantVector> const poly_vec){
        auto len = poly_vec.size();
        assert(not poly_vec.empty());
        if(len == 1){
            return poly_vec[0];
        }
        if(len == 2){
            return prod_sos(poly_vec[0],poly_vec[1]);
        }
        auto mid = len / 2;
        auto first_it = poly_vec.begin();
        auto mid_it= poly_vec.begin() + mid;
        auto last_it  = poly_vec.end();
        std::vector<InterpolantVector> vec1(first_it, mid_it);
        std::vector<InterpolantVector> vec2(mid_it, last_it);
        return prod_sos(prod_sos(vec1), prod_sos(vec2));
    }

    InterpolantVector generate_zero_polynomial();

    void add_polynomial(InterpolantVector &polynomial);

    Instance construct_SOS_instance();

    void print_solution(Solution sol);

    void plot_polynomials_and_solution(const Solution &sol);

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

};

#endif //NONSYMMETRICCONICOPTIMIZATION_ENVELOPEPROBLEMSOS_H
