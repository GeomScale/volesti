// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MAX_ELLIPSOID_ROUNDING_HPP
#define MAX_ELLIPSOID_ROUNDING_HPP

#include "max_inscribed_ellipsoid.hpp"
#include "analytic_center_linear_ineq.h"
#include "feasible_point.hpp"

enum EllipsoidType
{
  MAX_ELLIPSOID = 1,
  LOG_BARRIER = 2
};

template<int C>
struct inscribed_ellispoid
{
    template<typename MT, typename Custom_MT, typename VT, typename NT>
    inline static std::tuple<MT, VT, NT> 
    compute(Custom_MT A, VT b, VT const& x0,
            unsigned int const& maxiter,
            NT const& tol, NT const& reg) 
    {
      std::runtime_error("No roudning method is defined");
      return std::tuple<MT, VT, NT>();
    }
};

template <>
struct inscribed_ellispoid<EllipsoidType::MAX_ELLIPSOID>
{
    template<typename MT, typename Custom_MT, typename VT, typename NT>
    inline static std::tuple<MT, VT, NT> 
    compute(Custom_MT A, VT b, VT const& x0,
            unsigned int const& maxiter,
            NT const& tol, NT const& reg) 
    {
      return max_inscribed_ellipsoid<MT>(A, b, x0, maxiter, tol, reg);
    }
};

template <>
struct inscribed_ellispoid<EllipsoidType::LOG_BARRIER>
{
    template<typename MT, typename Custom_MT, typename VT, typename NT>
    inline static std::tuple<MT, VT, NT> 
    compute(Custom_MT const& A, VT const& b, VT const& x0,
            unsigned int const& maxiter,
            NT const& tol, NT&) 
    {
      return analytic_center_linear_ineq<MT, Custom_MT, VT, NT>(A, b, x0);
    }
};

template 
<
    typename MT,
    typename VT,
    typename NT,
    typename Polytope,
    int ellipsopid_type = EllipsoidType::MAX_ELLIPSOID
>
std::tuple<MT, VT, NT> max_inscribed_ellipsoid_rounding(Polytope &P, 
                                                        unsigned int const max_iterations = 5,
                                                        NT const max_eig_ratio = NT(6))
{
    typedef typename Polytope::PointType Point;
    VT x = compute_feasible_point(P.get_mat(), P.get_vec());
    return max_inscribed_ellipsoid_rounding<MT, VT, NT>(P, Point(x), max_iterations, max_eig_ratio);
}

template 
<
    typename MT,
    typename VT,
    typename NT,
    typename Polytope,
    typename Point,
    int ellipsopid_type = EllipsoidType::MAX_ELLIPSOID
>
std::tuple<MT, VT, NT> max_inscribed_ellipsoid_rounding(Polytope &P, 
                                                        Point const& InnerPoint,
                                                        unsigned int const max_iterations = 5,
                                                        NT const max_eig_ratio = NT(6))
{
    VT x0 = InnerPoint.getCoefficients(), center;
    MT E, L;
    bool converged;
    unsigned int maxiter = 500, iter = 1, d = P.dimension();

    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -4.0), round_val = 1.0;

    MT T = MT::Identity(d, d);
    VT shift = VT::Zero(d);

    while (true)
    {
        // compute the largest inscribed ellipsoid in P centered at x0
        //std::tie(E, center, converged) = max_inscribed_ellipsoid<MT>(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);
        std::tie(E, center, converged) = 
                       inscribed_ellispoid<ellipsopid_type>::template compute<MT>(P.get_mat(), P.get_vec(),
                                                                                  x0, maxiter, tol, reg);
        E = (E + E.transpose()) / 2.0;
        E += MT::Identity(d, d)*std::pow(10, -8.0); //normalize E

        Eigen::LLT<MT> lltOfA(E.llt().solve(MT::Identity(E.cols(), E.cols()))); // compute the Cholesky decomposition of E^{-1}
        L = lltOfA.matrixL();

        // computing eigenvalues of E
        Spectra::DenseSymMatProd<NT> op(E);
        // The value of ncv is chosen empirically
        Spectra::SymEigsSolver<NT, Spectra::SELECT_EIGENVALUE::BOTH_ENDS, 
                               Spectra::DenseSymMatProd<NT>> eigs(&op, 2, std::min(std::max(10, int(d)/5), int(d)));
        eigs.init();
        int nconv = eigs.compute();
        if (eigs.info() == Spectra::COMPUTATION_INFO::SUCCESSFUL) {
            R = 1.0 / eigs.eigenvalues().coeff(1);
            r = 1.0 / eigs.eigenvalues().coeff(0);
        } else  {
            Eigen::SelfAdjointEigenSolver<MT> eigensolver(E);
            if (eigensolver.info() == Eigen::ComputationInfo::Success) {
                R = 1.0 / eigensolver.eigenvalues().coeff(0);
                r = 1.0 / eigensolver.eigenvalues().template tail<1>().value();
            } else {
                std::runtime_error("Computations failed.");
            }
        }
        // shift polytope and apply the linear transformation on P
        P.shift(center);
        shift.noalias() += T * center;
        T.applyOnTheRight(L); // T = T * L;
        round_val *= L.transpose().determinant();
        P.linear_transformIt(L);

        reg = std::max(reg / 10.0, std::pow(10, -10.0));
        P.normalize();
        x0 = VT::Zero(d);

        // check the roundness of the polytope
        std::cout<<"std::abs(R / r): "<<std::abs(R / r)<<std::endl;
        if(((std::abs(R / r) <= max_eig_ratio && converged) || iter >= max_iterations)){
            break;
        }

        iter++;
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(round_val));
    return result;
}

#endif 
