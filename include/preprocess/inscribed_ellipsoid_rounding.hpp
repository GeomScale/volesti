// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef INSCRIBED_ELLIPSOID_ROUNDING_HPP
#define INSCRIBED_ELLIPSOID_ROUNDING_HPP

#include "preprocess/max_inscribed_ellipsoid.hpp"
#include "preprocess/barrier_center_ellipsoid.hpp"
#include "preprocess/feasible_point.hpp"


template<typename MT, int ellipsoid_type, typename Custom_MT, typename VT, typename NT>
inline static std::tuple<MT, VT, NT>
compute_inscribed_ellipsoid(Custom_MT A, VT b, VT const& x0,
                            unsigned int const& maxiter,
                            NT const& tol, NT const& reg)
{
    if constexpr (ellipsoid_type == EllipsoidType::MAX_ELLIPSOID)
    {
        return max_inscribed_ellipsoid<MT>(A, b, x0, maxiter, tol, reg);
    } else if constexpr (ellipsoid_type == EllipsoidType::LOG_BARRIER ||
                         ellipsoid_type == EllipsoidType::VOLUMETRIC_BARRIER ||
                         ellipsoid_type == EllipsoidType::VAIDYA_BARRIER)
    {
        return barrier_center_ellipsoid_linear_ineq<MT, ellipsoid_type, NT>(A, b, x0);
    } else
    {
        std::runtime_error("Unknown rounding method.");
    }
    return {};
}

template 
<
    typename MT,
    typename VT,
    typename NT,
    typename Polytope,
    int ellipsoid_type = EllipsoidType::MAX_ELLIPSOID
>
std::tuple<MT, VT, NT> inscribed_ellipsoid_rounding(Polytope &P, 
                                                    unsigned int const max_iterations = 5,
                                                    NT const max_eig_ratio = NT(6))
{
    typedef typename Polytope::PointType Point;
    VT x = compute_feasible_point(P.get_mat(), P.get_vec());
    return inscribed_ellipsoid_rounding<MT, VT, NT>(P, Point(x), max_iterations, max_eig_ratio);
}

template 
<
    typename MT,
    typename VT,
    typename NT,
    typename Polytope,
    typename Point,
    int ellipsoid_type = EllipsoidType::MAX_ELLIPSOID
>
std::tuple<MT, VT, NT> inscribed_ellipsoid_rounding(Polytope &P, 
                                                    Point const& InnerPoint,
                                                    unsigned int const max_iterations = 5,
                                                    NT const max_eig_ratio = NT(6))
{
    unsigned int maxiter = 500, iter = 1, d = P.dimension();
    VT x0 = InnerPoint.getCoefficients(), center, shift = VT::Zero(d);
    MT E, L, T = MT::Identity(d, d);
    bool converged;
    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -4.0), round_val = 1.0;

    while (true)
    {
        // Compute the desired inscribed ellipsoid in P
        std::tie(E, center, converged) = 
            compute_inscribed_ellipsoid<MT, ellipsoid_type>(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);
        
        E = (E + E.transpose()) / 2.0;
        E += MT::Identity(d, d)*std::pow(10, -8.0); //normalize E

        Eigen::LLT<MT> lltOfA(E.llt().solve(MT::Identity(E.cols(), E.cols()))); // compute the Cholesky decomposition of E^{-1}
        L = lltOfA.matrixL();

        // Computing eigenvalues of E
        Spectra::DenseSymMatProd<NT> op(E);
        // The value of ncv is chosen empirically
        Spectra::SymEigsSolver<NT, Spectra::SELECT_EIGENVALUE::BOTH_ENDS, 
                               Spectra::DenseSymMatProd<NT>> eigs(&op, 2, std::min(std::max(10, int(d)/5), int(d)));
        eigs.init();
        int nconv = eigs.compute();
        if (eigs.info() == Spectra::COMPUTATION_INFO::SUCCESSFUL) {
            R = 1.0 / eigs.eigenvalues().coeff(1);
            r = 1.0 / eigs.eigenvalues().coeff(0);
        } else {
            Eigen::SelfAdjointEigenSolver<MT> eigensolver(E);
            if (eigensolver.info() == Eigen::ComputationInfo::Success) {
                R = 1.0 / eigensolver.eigenvalues().coeff(0);
                r = 1.0 / eigensolver.eigenvalues().template tail<1>().value();
            } else {
                std::runtime_error("Computations failed.");
            }
        }
        // Shift polytope and apply the linear transformation on P
        P.shift(center);
        shift.noalias() += T * center;
        T.applyOnTheRight(L); // T = T * L;
        round_val *= L.transpose().determinant();
        P.linear_transformIt(L);

        reg = std::max(reg / 10.0, std::pow(10, -10.0));
        P.normalize();
        x0 = VT::Zero(d);

        // Check the roundness of the polytope
        if(((std::abs(R / r) <= max_eig_ratio && converged) || iter >= max_iterations)) {
            break;
        }

        iter++;
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(round_val));
    return result;
}

#endif 
