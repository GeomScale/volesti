#ifndef SAMPLERS_RANDOM_POINT_GEN_MULTI_HPP
#define SAMPLERS_RANDOM_POINT_GEN_MULTI_HPP

#include "diagnostics/effective_sample_size.hpp"

template
<
    typename Walk
>
struct RandomPointEfficientGenerator
{
    template
    <
        typename Polytope,
        typename VT,
        typename MT,
        typename RandomNumberGenerator,
        typename Parameters
    >
    static void apply(Polytope& P,
                      VT &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      unsigned int &window,
                      unsigned int &Neff_sampled,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool &complete,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        typedef double NT;
        Walk walk(P, p, rng, parameters);
        ESSestimator<NT, VT, MT> estimator(window, P.dimension());
        //int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        walk.template burnin(P, p, nburns, walk_length, rng);
        bool done = false;
        unsigned int pointer = 0, total_ef_samples = 0, total_samples = 0, points_to_sample = rnum;
        int min_eff_samples, min_window = window, iter_points = 0;
        MT winPoints(P.dimension(), window);
        
        while (!done) 
        {
            for (int i = 0; i < min_window; i++)
            {
                walk.template apply(P, walk_length, rng);
                winPoints.col(pointer + i) = walk.template get_curr_sample();
                //if ((i+1)%100 == 0) {
                    //std::cout<<"number of sample points = "<<i<<std::endl;
                //}
            }
            estimator.update_estimator(winPoints);
            total_samples += min_window;
            if (total_samples >= TotalRandPoints.cols()) {
                if (total_samples > TotalRandPoints.cols()) {
                    TotalRandPoints.conservativeResize(total_samples, P.dimension());
                }
                done = true;
            }
            TotalRandPoints.block(total_samples - min_window, 0, min_window, P.dimension()) = winPoints.transpose();
            //min_eff_samples = int(ess_univariate_fft<NT, VT>(winPoints).minCoeff());
            //std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (done || total_samples => points_to_sample) {
                estimator.estimate_effective_sample_size();                       
                min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
                if (done && min_eff_samples < rnum) {
                    Neff_sampled = min_eff_samples;
                    return;
                }
                if (min_eff_samples >= rnum) {
                    complete = true;
                    Neff_sampled = min_eff_samples
                    return;
                }
                points_to_sample = (total_samples / min_eff_samples) * (rnum - min_eff_samples) + 200;
            }
        }
    }

    template
    <
            typename Polytope,
            typename VT,
            typename MT,
            typename RandomNumberGenerator
    >
    static void apply(Polytope& P,
                      VT &p,   // a point to start
                      unsigned int const& rnum,
                      unsigned int const& walk_length,
                      unsigned int &window,
                      unsigned int &Neff_sampled,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool &complete,
                      RandomNumberGenerator &rng)
    {
        typedef double NT;
        Walk walk(P, p, rng);
        ESSestimator<NT, VT, MT> estimator(window, P.dimension());
        //int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        walk.template burnin(P, p, nburns, walk_length, rng);
        bool done = false;
        unsigned int pointer = 0, total_ef_samples = 0, total_samples = 0, points_to_sample = rnum;
        int min_eff_samples, min_window = window, iter_points = 0;
        MT winPoints(P.dimension(), window);
        
        while (!done) 
        {
            for (int i = 0; i < min_window; i++)
            {
                walk.template apply(P, walk_length, rng);
                winPoints.col(pointer + i) = walk.template get_curr_sample();
                //if ((i+1)%100 == 0) {
                    //std::cout<<"number of sample points = "<<i<<std::endl;
                //}
            }
            estimator.update_estimator(winPoints);
            total_samples += min_window;
            if (total_samples >= TotalRandPoints.cols()) {
                if (total_samples > TotalRandPoints.cols()) {
                    TotalRandPoints.conservativeResize(total_samples, P.dimension());
                }
                done = true;
            }
            TotalRandPoints.block(total_samples - min_window, 0, min_window, P.dimension()) = winPoints.transpose();
            //min_eff_samples = int(ess_univariate_fft<NT, VT>(winPoints).minCoeff());
            //std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (done || total_samples => points_to_sample) {
                estimator.estimate_effective_sample_size();                       
                min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
                if (done && min_eff_samples < rnum) {
                    Neff_sampled = min_eff_samples;
                    return;
                }
                if (min_eff_samples >= rnum) {
                    complete = true;
                    Neff_sampled = min_eff_samples
                    return;
                }
                points_to_sample = (total_samples / min_eff_samples) * (rnum - min_eff_samples) + 200;
            }
        }
    }
};

#endif
