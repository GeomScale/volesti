// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef MAX_ELLIPSOID_ROUNDING_HPP
#define MAX_ELLIPSOID_ROUNDING_HPP

#include "max_inscribed_ellipsoid.hpp"

template 
<
    typename MT,
    typename VT,
    typename NT,
    typename Polytope,
    typename Point   
>
std::tuple<MT, VT, NT> max_inscribed_ellipsoid_rounding(Polytope &P, 
                                                        Point const& InnerPoint,
                                                        unsigned int const max_iterations = 5,
                                                        NT const max_eig_ratio = NT(6))
{
    std::pair<std::pair<MT, VT>, bool> iter_res;
    iter_res.second = false;

    VT x0 = InnerPoint.getCoefficients();
    MT E, L;
    unsigned int maxiter = 500, iter = 1, d = P.dimension();

    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -4.0), round_val = 1.0;

    MT T = MT::Identity(d, d);
    VT shift = VT::Zero(d);

    while (true)
    {
        // compute the largest inscribed ellipsoid in P centered at x0
        iter_res = max_inscribed_ellipsoid<MT>(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);
        E = iter_res.first.first;
        E = (E + E.transpose()) / 2.0;
        E = E + MT::Identity(d, d)*std::pow(10, -8.0); //normalize E

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
        P.shift(iter_res.first.second);
        shift += T * iter_res.first.second;
        T = T * L;
        round_val *= L.transpose().determinant();
        P.linear_transformIt(L);

        reg = std::max(reg / 10.0, std::pow(10, -10.0));
        P.normalize();
        x0 = VT::Zero(d);
        
        // check the roundness of the polytope
        if(((std::abs(R / r) <= max_eig_ratio && iter_res.second) || iter >= max_iterations)){
            break;
        }

        iter++;
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(round_val));
    return result;
}

#endif 
