// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_ORACLE_AUTODIFF_FUNCTORS_HPP
#define ODE_SOLVERS_ORACLE_AUTODIFF_FUNCTORS_HPP

#include "Eigen/Eigen"
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include "cartesian_geom/cartesian_kernel.h"
#include "cartesian_geom/autopoint.h"

struct AutoDiffFunctor {

    template <typename NT>
    struct parameters {
        unsigned int order;
        NT L; // Lipschitz constant for gradient
        NT m; // Strong convexity constant
        NT kappa; // Condition number
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> data;
        parameters() : order(2), L(4), m(4), kappa(1){};
    };

    template <typename NT>
    struct FunctionFunctor_internal {
        using Autopoint = autopoint<NT>;
        using Coeff = typename autopoint<NT>::Coeff;
        using FT = typename autopoint<NT>::FT;
        using Point = typename Cartesian<NT>::Point;

        static std::function<FT(const Autopoint &, const Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> &)> pdf;

        FT static result_internal(const Coeff &x, const Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> &data){
            return pdf(x, data); //
        }

        // external interface
        Point static differentiate(Point const &x0, const Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> &data) {
            Autopoint x = Autopoint(x0.getCoefficients()); // cast into autopoint
            auto x1 = x.getCoefficients();
            Coeff y = autodiff::gradient(result_internal, autodiff::wrt(x1), autodiff::at(x1, data));
            auto result = y.template cast<NT>();
            return -1 * Point(result);
        }

        NT static result(Point const &x0, const Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> &data) {
            Autopoint x = Autopoint(x0.getCoefficients()); // cast to autopoint
            auto x1 = x.getCoefficients();
            return result_internal(x1, data).val();
        }
    };

    template <typename Point>
    struct GradientFunctor {
        using NT = typename Point::FT;

        FunctionFunctor_internal<NT> F;
        parameters<NT> &params;

        GradientFunctor(parameters<NT> &params_) : params(params_){};

        // The index i represents the state vector index
        Point operator()(unsigned int const &i, std::vector<Point> const &xs, NT const &t) const {
            // std::cout<<"calling gradient functor"<<std::flush;
            if (i == params.order - 1) {
                return F.differentiate(xs[0], params.data);
            }
            else {
                return xs[i + 1]; // returns derivative
            }
        }
    };

    template <typename Point>
    struct FunctionFunctor {
        using NT = typename Point::FT;
        parameters<NT> &params;

        FunctionFunctor(parameters<NT> &params_) : params(params_){};
        // The index i represents the state vector index
        FunctionFunctor_internal<NT> F;

        NT operator()(Point const &x) const {
            return F.result(x, params.data);
        }
    };
};

#endif //ODE_SOLVERS_ORACLE_AUTODIFF_FUNCTORS_HPP