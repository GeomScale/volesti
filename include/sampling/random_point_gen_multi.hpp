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
                      unsigned int const& max_num_samples,
                      unsigned int const& walk_length,
                      unsigned int &window,
                      unsigned int &Neff_sampled,
                      unsigned int &total_samples,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool &complete,
                      RandomNumberGenerator &rng,
                      bool request_rounding,
                      Parameters const& parameters)
                      //total_num_sampled TODO!!!
    {
        typedef double NT;
        Walk walk(P, p, rng, parameters);
        ESSestimator<NT, VT, MT> estimator(window, P.dimension());
        //int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        NT avg_ref = 0.0;
        //MT starting_points(P.dimension(), 10+int(std::log(NT(P.dimension()))));
        walk.template burnin(P, p, 10+int(std::log(NT(P.dimension()))), 10, avg_ref, rng);
        std::cout<<"avg_ref = "<<avg_ref<<std::endl;
        bool done = false;
        unsigned int points_to_sample = rnum;
        int min_eff_samples, min_window = window;
        total_samples = 0;
        MT winPoints(P.dimension(), window);
        std::cout<<"request_rounding = "<<request_rounding<<std::endl;
        std::cout<<"rnum = "<<rnum<<std::endl;
        VT q(P.dimension());
        while (!done) 
        {
            std::cout<<"[j] points_to_sample = "<<points_to_sample<<std::endl;
            //std::cout<<"j = "<<j<<std::endl;
            walk.template get_starting_point(P, p, q, 10, rng);
            walk.template initialize(P, q, rng);
            for (int i = 0; i < min_window; i++)
            {
                walk.template apply(P, walk_length, rng);
                winPoints.col(i) = walk.template get_curr_sample();
                if ((i+1)%100 == 0) {
                    std::cout<<"number of sample points = "<<i<<std::endl;
                }
            }
            estimator.update_estimator(winPoints);
            total_samples += min_window;
            std::cout<<"total_samples = "<<total_samples<<std::endl;
            if (total_samples >= TotalRandPoints.rows()) {
                if (total_samples > TotalRandPoints.rows()) {
                    TotalRandPoints.conservativeResize(total_samples, P.dimension());
                }
                if (request_rounding || total_samples >= max_num_samples) done = true;
            }
            std::cout<<"done = "<<done<<std::endl;
            TotalRandPoints.block(total_samples - min_window, 0, min_window, P.dimension()) = winPoints.transpose();
            //min_eff_samples = int(ess_univariate_fft<NT, VT>(winPoints).minCoeff());
            //std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (done || total_samples >= points_to_sample) {
                estimator.estimate_effective_sample_size();                       
                min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
                std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
                if (done && min_eff_samples < rnum) {
                    Neff_sampled = min_eff_samples;
                    return;
                }
                if (min_eff_samples >= rnum) {
                    complete = true;
                    std::cout<<"complete = "<<complete<<std::endl;
                    Neff_sampled = min_eff_samples;
                    return;
                }
                if (min_eff_samples > 0) {
                    points_to_sample += (total_samples / min_eff_samples) * (rnum - min_eff_samples) + 100;
                } else {
                    points_to_sample = total_samples + 100;
                }
                std::cout<<"points_to_sample = "<<points_to_sample<<std::endl;
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
                      unsigned int const& max_num_samples,
                      unsigned int const& walk_length,
                      unsigned int &window,
                      unsigned int &Neff_sampled,
                      unsigned int &total_samples,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool &complete,
                      bool request_rounding,
                      RandomNumberGenerator &rng)
    {
        typedef double NT;
        Walk walk(P, p, rng);
        ESSestimator<NT, VT, MT> estimator(window, P.dimension());
        //int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        NT avg_ref = 0.0;
        //MT starting_points(P.dimension(), 10+int(std::log(NT(P.dimension()))));
        walk.template burnin(P, p, 10+int(std::log(NT(P.dimension()))), 10, avg_ref, rng);
        std::cout<<"avg_ref = "<<avg_ref<<std::endl;
        bool done = false;
        unsigned int points_to_sample = rnum;
        int min_eff_samples, min_window = window;
        total_samples = 0;
        MT winPoints(P.dimension(), window);
        std::cout<<"request_rounding = "<<request_rounding<<std::endl;
        std::cout<<"rnum = "<<rnum<<std::endl;
        VT q(P.dimension());
        while (!done) 
        {
            std::cout<<"[j] points_to_sample = "<<points_to_sample<<std::endl;
            //std::cout<<"j = "<<j<<std::endl;
            walk.template get_starting_point(P, p, q, 10, rng);
            walk.template initialize(P, q, rng);
            for (int i = 0; i < min_window; i++)
            {
                walk.template apply(P, walk_length, rng);
                winPoints.col(i) = walk.template get_curr_sample();
                if ((i+1)%100 == 0) {
                    std::cout<<"number of sample points = "<<i<<std::endl;
                }
            }
            estimator.update_estimator(winPoints);
            total_samples += min_window;
            std::cout<<"total_samples = "<<total_samples<<std::endl;
            if (total_samples >= TotalRandPoints.rows()) {
                if (total_samples > TotalRandPoints.rows()) {
                    TotalRandPoints.conservativeResize(total_samples, P.dimension());
                }
                if (request_rounding || total_samples >= max_num_samples) done = true;
            }
            std::cout<<"done = "<<done<<std::endl;
            TotalRandPoints.block(total_samples - min_window, 0, min_window, P.dimension()) = winPoints.transpose();
            //min_eff_samples = int(ess_univariate_fft<NT, VT>(winPoints).minCoeff());
            //std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (done || total_samples >= points_to_sample) {
                estimator.estimate_effective_sample_size();                       
                min_eff_samples = int(estimator.get_effective_sample_size().minCoeff());
                std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
                if (done && min_eff_samples < rnum) {
                    Neff_sampled = min_eff_samples;
                    return;
                }
                if (min_eff_samples >= rnum) {
                    complete = true;
                    std::cout<<"complete = "<<complete<<std::endl;
                    Neff_sampled = min_eff_samples;
                    return;
                }
                if (min_eff_samples > 0) {
                    points_to_sample += (total_samples / min_eff_samples) * (rnum - min_eff_samples) + 100;
                } else {
                    points_to_sample = total_samples + 100;
                }
                std::cout<<"points_to_sample = "<<points_to_sample<<std::endl;
            }
        }
    }
};

#endif
