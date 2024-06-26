// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef SVD_ROUNDING_HPP
#define SVD_ROUNDING_HPP

#include "diagnostics/univariate_psrf.hpp"
#include "sampling/sampling.hpp"
#include "max_inscribed_ellipsoid.hpp"

template
<
    typename WalkTypePolicy,
    int simdLen = 1,
    typename Polytope,
    typename Point,
    typename MT,
    typename VT,
    typename RandomNumberGenerator
>
void svd_on_sample(Polytope &P, Point &p, unsigned int const& num_rounding_steps, MT &V, VT &s, VT &Means,
                   unsigned int const& walk_length, RandomNumberGenerator &rng)
{


    using NT = double;

    std::pair<std::pair<MT, VT>, bool> iter_res;
    iter_res.second = false;
    
    MT E, L;
    unsigned int maxiter = 500, iter = 1, d = P.dimension();
    const Point X = p;
    VT x0 = X.getCoefficients();

    NT R = 100.0, r = 1.0, tol = std::pow(10, -6.0), reg = std::pow(10, -4.0), round_val = 1.0;

    MT T = MT::Identity(d, d);
    VT shift = VT::Zero(d);

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

    std::cout << "Current ratio of max inscribed ellipsoid R/r " << R / r << std::endl;








    std::list<Point> randPoints;

    unsigned int N = num_rounding_steps;

    std::cout << "Sampling " << N << " points" << std::endl;

    //P.print();

    if constexpr (std::is_same_v<WalkTypePolicy, CRHMCWalk>) {

        std::pair<Point, NT> InnerBall = P.ComputeInnerBall();

        using Func = GaussianFunctor::FunctionFunctor<Point>;
        using Grad = GaussianFunctor::GradientFunctor<Point>;
        using Hess = GaussianFunctor::HessianFunctor<Point>;
        using func_params=GaussianFunctor::parameters<NT, Point>;
        func_params params = func_params(InnerBall.first, 1e-6, 1);
        Func* f= new Func(params);
        Grad* g= new Grad(params);
        Hess* h= new Hess(params);
        
        unsigned int nr_burns = N;

        execute_crhmc<Polytope, RandomNumberGenerator, std::list<Point>, Grad, Func, Hess, WalkTypePolicy, simdLen>(
        P, rng, randPoints, walk_length, N, nr_burns, g, f, h);

        delete f;
        delete g;
        delete h;







        bool ok = true;

        for(auto pp:randPoints)
        {
            if(P.is_in(pp) != -1)
                ok = false;
            if(pp.dimension() != P.dimension())
                ok = false;
        }

        MT samples = MT(randPoints.front().dimension(), randPoints.size());
        int i=0;
        for (typename std::list<Point>::iterator it = randPoints.begin(); it != randPoints.end(); ++it){
            samples.col(i) = (*it).getCoefficients();
            i++;
        }
        MT max_psrf = univariate_psrf<NT, VT>(samples);

        std::cout<<"PSRF: "<<max_psrf<<std::endl;
        std::cout << ok << std::endl;








    } else {

        typedef typename WalkTypePolicy::template Walk
                <
                        Polytope,
                        RandomNumberGenerator
                > walk;

        typedef RandomPointGenerator <walk> RandomPointGenerator;
        PushBackWalkPolicy push_back_policy;

        RandomPointGenerator::apply(P, p, N, walk_length, randPoints,
                                    push_back_policy, rng);
    }

    MT RetMat(N, P.dimension());

    int jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        RetMat.row(jj) = (*rpit).getCoefficients().transpose();
    }

    for (int i = 0; i < P.dimension(); ++i) {
        Means(i) = RetMat.col(i).mean();
    }

    for (int i = 0; i < N; ++i) {
        RetMat.row(i) = RetMat.row(i) - Means.transpose();
    }

    Eigen::BDCSVD<MT> svd(RetMat, Eigen::ComputeFullV);
    s = svd.singularValues() / svd.singularValues().minCoeff();

    if (s.maxCoeff() >= 2.0) {
        for (int i = 0; i < s.size(); ++i) {
            if (s(i) < 2.0) {
                s(i) = 1.0;
            }
        }
        V = svd.matrixV();
    } else {
        s = VT::Ones(P.dimension());
        V = MT::Identity(P.dimension(), P.dimension());
    }
}


template
<
    typename WalkTypePolicy,
    typename MT,
    typename VT,
    int simdLen = 1,
    typename Polytope,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>

std::tuple<MT, VT, NT> svd_rounding(Polytope &P,
                                    std::pair<Point,NT> &InnerBall,
                                    const unsigned int &walk_length,
                                    RandomNumberGenerator &rng)
{
    NT tol = 0.00000001;
    NT R = std::pow(10,10), r = InnerBall.second;

    int n = P.dimension(), m = P.num_of_hyperplanes();

    Polytope P_old(P);

    MT old_A(m,n), A = P.get_mat(), T = MT::Identity(n,n), round_mat, r_inv;
    VT T_shift = VT::Zero(n), shift(n), s(n);

    Point p(n);

    bool done = false, last_round_under_p, fail;

    unsigned int tries=0, num_rounding_steps = 10 * n, rounding_samples = 0, round_it;
    NT max_s, s_cutof, p_cutof, num_its, prev_max_s = std::numeric_limits<NT>::max(),
       s_cutoff, p_cutoff;
    MT V(n,n), S(n,n);

    while (!done) {

        T = MT::Identity(n, n);
        T_shift = VT::Zero(n);

        round_it = 1;
        max_s = std::numeric_limits<NT>::max();
        s_cutoff = 2.3;
        p_cutoff = 10.0;
        last_round_under_p = false;
        fail = false;
        num_its = 20;

        while (max_s > s_cutoff && round_it <= num_its) {

            p = InnerBall.first;
            svd_on_sample<WalkTypePolicy, simdLen>(P, p, num_rounding_steps, V, s,
                                          shift, walk_length, rng);

            rounding_samples = rounding_samples + num_rounding_steps;
            max_s = s.maxCoeff();

            if (max_s <= p_cutoff && max_s > s_cutoff) {
                if (last_round_under_p) {
                    num_rounding_steps = num_rounding_steps * 2;
                    p = InnerBall.first;
                    svd_on_sample<WalkTypePolicy, simdLen>(P, p, num_rounding_steps, V, s,
                                                  shift, walk_length, rng);
                    max_s = s.maxCoeff();
                } else {
                    last_round_under_p = true;
                }
            } else {
                last_round_under_p = false;
            }
            S = s.asDiagonal();
            round_mat = V * S;
            r_inv = VT::Ones(n).cwiseProduct(s.cwiseInverse()).asDiagonal() * V.transpose();

            if (round_it != 1 && max_s >= NT(4) * prev_max_s) {
                fail = true;
                break;
            }

            round_it++;
            prev_max_s = max_s;

            P.shift(shift);
            P.linear_transformIt(round_mat);
            InnerBall = P.ComputeInnerBall();
            T_shift += T * shift;
            T = T * round_mat;
        }
        if (round_it <= num_its && !fail) {
            done = true;
        } else {
            tries = tries + 1;
            num_rounding_steps = num_rounding_steps * 2.0;
            P = P_old;
        }
    }

    std::tuple<MT, VT, NT> result = std::make_tuple(T, shift, std::abs(T.determinant()));
    return result;
}


#endif
