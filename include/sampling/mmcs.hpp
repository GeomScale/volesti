// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MMCS_HPP
#define MMCS_HPP

#include "diagnostics/ess_window_updater.hpp"

/**
 *  The class implements a single step of the Multiphase Monte Carlo Sampling algorithm
 *  given in,
 *
 *  A. Chalkis, V. Fisikopoulos, E. Tsigaridas, H. Zafeiropoulos, Geometric algorithms for sampling the flux space of metabolic networks, SoCG 21.
 *
 * @tparam Polytope convex polytope type
 * @tparam RandomNumberGenerator random number generator type
 * @tparam MT matrix type
 * @tparam Point cartensian point type
 * @tparam WalkTypePolicy random walk type
*/
template
<
    typename Polytope,
    typename RandomNumberGenerator,
    typename MT,
    typename Point,
    typename WalkTypePolicy
>
bool perform_mmcs_step(Polytope &P,
                       RandomNumberGenerator &rng,
                       unsigned int const& walk_length,
                       unsigned int const& target_ess,
                       unsigned int const& max_num_samples,
                       unsigned int const& window,
                       unsigned int &Neff_sampled,
                       unsigned int &total_samples,
                       unsigned int const& num_rounding_steps,
                       MT &TotalRandPoints,
                       const Point &starting_point,
                       unsigned int const& nburns,
                       bool request_rounding,
                       WalkTypePolicy &WalkType)
{
    typedef typename Polytope::NT NT;
    typedef typename Polytope::VT VT;
    typedef typename WalkTypePolicy::template Walk
    <
        Polytope,
        RandomNumberGenerator
    > Walk;

    bool done = false;
    unsigned int points_to_sample = target_ess;
    int min_eff_samples;
    total_samples = 0;
    MT winPoints(P.dimension(), window);
    Point q(P.dimension());

    Point p = starting_point;

    if (request_rounding)
    {
        TotalRandPoints.setZero(num_rounding_steps, P.dimension());
    }
    else
    {
        TotalRandPoints.setZero(max_num_samples, P.dimension());
    }

    Walk walk(P, p, rng, WalkType.param);
    ESSestimator<NT, VT, MT> estimator(window, P.dimension());

    walk.template parameters_burnin(P, p, 10 + int(std::log(NT(P.dimension()))), 10, rng);

    while (!done)
    {
        walk.template get_starting_point(P, p, q, 10, rng);
        for (int i = 0; i < window; i++)
        {
            walk.apply(P, q, walk_length, rng);
            winPoints.col(i) = q.getCoefficients();
        }
        estimator.update_estimator(winPoints);
        total_samples += window;
        if (total_samples >= TotalRandPoints.rows())
        {
            if (total_samples > TotalRandPoints.rows())
            {
                TotalRandPoints.conservativeResize(total_samples, P.dimension());
            }
            if (request_rounding || total_samples >= max_num_samples)
            {
                done = true;
            }
        }
        TotalRandPoints.block(total_samples - window, 0, window, P.dimension()) = winPoints.transpose();
        if (done || total_samples >= points_to_sample)
        {
            estimator.estimate_effective_sample_size();
            min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
            if (done && min_eff_samples < target_ess)
            {
                Neff_sampled = min_eff_samples;
                return false;
            }
            if (min_eff_samples >= target_ess)
            {
                Neff_sampled = min_eff_samples;
                return true;
            }
            if (min_eff_samples > 0)
            {
                points_to_sample += (total_samples / min_eff_samples) * (target_ess - min_eff_samples) + 100;
            }
            else
            {
                points_to_sample = total_samples + 100;
            }
        }
    }
    return false;
}

template
<
    typename Polytope,
    typename MT
>
void mmcs(Polytope const& Pin,
          int const& Neff,
          MT& S,
          int& total_neff)
{
    mmcs(Pin, Neff, S, total_neff, 1);
}

template
<
    typename Polytope,
    typename MT
>
void mmcs(Polytope const& Pin,
          int const& Neff,
          MT& S,
          int& total_neff,
          unsigned int const& walk_length)
{
    using RNGType = BoostRandomNumberGenerator<boost::mt19937, typename Polytope::NT>;
    RNGType rng(Pin.dimension());
    mmcs(Pin, Neff, S, total_neff, walk_length, rng);
}

template
<
    typename Polytope,
    typename MT,
    typename RandomNumberGenerator
>
void mmcs(Polytope const& Pin,
          int const& Neff,
          MT& S,
          int& total_neff,
          unsigned int const& walk_length,
          RandomNumberGenerator &rng)
{
    using NT = typename Polytope::NT;
    using VT = typename Polytope::VT;
    using Point = typename Polytope::PointType;

    auto P = Pin;
    const int n = P.dimension();
    MT T = MT::Identity(n, n);
    VT T_shift = VT::Zero(n);

    unsigned int current_Neff = Neff;
    unsigned int round_it = 1;
    unsigned int num_rounding_steps = 20 * n;
    unsigned int num_its = 20;
    unsigned int phase = 0;
    unsigned int window = 100;
    unsigned int max_num_samples = 100 * n;
    unsigned int total_samples;
    unsigned int nburns;
    unsigned int total_number_of_samples_in_P0 = 0;

    NT max_s;
    NT s_cutoff = 3.0;
    NT L;
    bool complete = false;
    bool request_rounding = true;
    bool rounding_completed = false;
    bool req_round_temp = request_rounding;

    std::pair<Point, NT> InnerBall;

    std::cout << "target effective sample size = " << current_Neff << "\n" << std::endl;

    while(true)
    {
        phase++;
        if (request_rounding && rounding_completed)
        {
            req_round_temp = false;
        }

        if (req_round_temp)
        {
            nburns = num_rounding_steps / window + 1;
        }
        else
        {
            nburns = max_num_samples / window + 1;
        }

        InnerBall = P.ComputeInnerBall();
        L = NT(6) * std::sqrt(NT(n)) * InnerBall.second;
        AcceleratedBilliardWalk WalkType(L);

        unsigned int Neff_sampled;
        MT TotalRandPoints;
        complete = perform_mmcs_step(P, rng, walk_length, current_Neff, max_num_samples, window,
                                     Neff_sampled, total_samples, num_rounding_steps, TotalRandPoints,
                                     InnerBall.first, nburns, req_round_temp, WalkType);

        current_Neff -= Neff_sampled;
        std::cout << "phase " << phase << ": number of correlated samples = " << total_samples << ", effective sample size = " << Neff_sampled;
        total_neff += Neff_sampled;
        Neff_sampled = 0;

        MT Samples = TotalRandPoints.transpose(); //do not copy TODO!
        for (int i = 0; i < total_samples; i++)
        {
            Samples.col(i) = T * Samples.col(i) + T_shift;
        }

        S.conservativeResize(P.dimension(), total_number_of_samples_in_P0 + total_samples);
        S.block(0, total_number_of_samples_in_P0, P.dimension(), total_samples) = Samples.block(0, 0, P.dimension(), total_samples);
        total_number_of_samples_in_P0 += total_samples;
        if (!complete)
        {
            if (request_rounding && !rounding_completed)
            {
                VT shift(n), s(n);
                MT V(n,n), S(n,n), round_mat;
                for (int i = 0; i < P.dimension(); ++i)
                {
                    shift(i) = TotalRandPoints.col(i).mean();
                }

                for (int i = 0; i < total_samples; ++i)
                {
                    TotalRandPoints.row(i) = TotalRandPoints.row(i) - shift.transpose();
                }

                Eigen::BDCSVD<MT> svd(TotalRandPoints, Eigen::ComputeFullV);
                s = svd.singularValues() / svd.singularValues().minCoeff();

                if (s.maxCoeff() >= 2.0)
                {
                    for (int i = 0; i < s.size(); ++i)
                    {
                        if (s(i) < 2.0)
                        {
                            s(i) = 1.0;
                        }
                    }
                    V = svd.matrixV();
                }
                else
                {
                    s = VT::Ones(P.dimension());
                    V = MT::Identity(P.dimension(), P.dimension());
                }
                max_s = s.maxCoeff();
                S = s.asDiagonal();
                round_mat = V * S;

                round_it++;
                P.shift(shift);
                P.linear_transformIt(round_mat);
                T_shift += T * shift;
                T = T * round_mat;

                std::cout << ", ratio of the maximum singilar value over the minimum singular value = " << max_s << std::endl;

                if (max_s <= s_cutoff || round_it > num_its)
                {
                    rounding_completed = true;
                }
            }
            else
            {
                std::cout<<"\n";
            }
        }
        else
        {
            std::cout<<"\n\n";
            break;
        }
    }
}

#endif
