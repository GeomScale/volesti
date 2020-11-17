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
                      unsigned int &min_skip,
                      unsigned int &window,
                      MT &EssRandPoints,
                      unsigned int &Neff_sampled,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool rounding_requested,
                      bool &complete,
                      RandomNumberGenerator &rng,
                      Parameters const& parameters)
    {
        typedef double NT;
        Walk walk(P, p, rng, parameters);
        //int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        walk.template burnin(P, p, nburns, walk_length, rng);
        bool done = false;
        unsigned int pointer = 0, total_ef_samples = 0, total_samples = 0;
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
            total_samples += min_window;
            min_eff_samples = int(ess_univariate_fft<NT, VT>(winPoints).minCoeff());
            //std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (min_eff_samples <= 0) {
                std::cout<<"[Complete] min_eff_samples = "<<min_eff_samples<<std::endl;
                pointer = winPoints.cols();
                winPoints.conservativeResize(P.dimension(), winPoints.cols() + window);
                std::cout<<"window  = "<<winPoints.cols()<<std::endl;
                min_window = window;
                //return;
            } else {
                total_ef_samples += min_eff_samples;
                //std::cout<<"total_ef_samples  = "<<total_ef_samples<<std::endl;
                std::cout<<"winPoints.cols()  = "<<winPoints.cols()<<std::endl;
                if (winPoints.cols() / min_eff_samples > min_skip) {
                    min_skip = winPoints.cols() / min_eff_samples;
                }
                std::cout<<"min_skip  = "<<min_skip<<std::endl;
                for (int i = min_skip-2; i < winPoints.cols(); i+=min_skip)
                {
                    //std::cout<<"indices in window added  = "<<i<<std::endl;
                    EssRandPoints.col(iter_points) = winPoints.col(i);
                    iter_points++;
                    if (iter_points == rnum) {
                        Neff_sampled = iter_points;
                        complete = true;
                        return;
                    }
                }
                Neff_sampled = iter_points;
                pointer = 0;
                if (rounding_requested) {
                    if (total_samples > TotalRandPoints.cols()) {
                        TotalRandPoints.conservativeResize(total_samples, P.dimension());
                        done = true;
                    }
                    TotalRandPoints.block(total_samples - min_window, 0, min_window, P.dimension()) = winPoints.transpose();
                }
                min_window = winPoints.cols();                
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
                      unsigned int &min_skip,
                      unsigned int &window,
                      MT &EssRandPoints,
                      unsigned int &Neff_sampled,
                      MT &TotalRandPoints,
                      unsigned int const& nburns,
                      bool rounding_requested,
                      bool &complete,
                      RandomNumberGenerator &rng)
    {
        typedef double NT;
        Walk walk(P, p, rng);
        //int num_points = 100 + 2*int( std::sqrt(NT(P.dimension())) );
        walk.template burnin(P, p, nburns, walk_length, rng);
        bool done = false;
        unsigned int pointer = 0, total_ef_samples = 0, total_samples = 0;
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
            total_samples += min_window;
            min_eff_samples = int(ess_univariate_fft<NT, VT>(winPoints).minCoeff());
            //std::cout<<"min_eff_samples = "<<min_eff_samples<<std::endl;
            if (min_eff_samples <= 0) {
                std::cout<<"[Complete] min_eff_samples = "<<min_eff_samples<<std::endl;
                pointer = winPoints.cols();
                winPoints.conservativeResize(P.dimension(), winPoints.cols() + window);
                std::cout<<"window  = "<<winPoints.cols()<<std::endl;
                min_window = window;
                //return;
            } else {
                total_ef_samples += min_eff_samples;
                //std::cout<<"total_ef_samples  = "<<total_ef_samples<<std::endl;
                std::cout<<"winPoints.cols()  = "<<winPoints.cols()<<std::endl;
                if (winPoints.cols() / min_eff_samples > min_skip) {
                    min_skip = winPoints.cols() / min_eff_samples;
                }
                std::cout<<"min_skip  = "<<min_skip<<std::endl;
                for (int i = min_skip-2; i < winPoints.cols(); i+=min_skip)
                {
                    //std::cout<<"indices in window added  = "<<i<<std::endl;
                    EssRandPoints.col(iter_points) = winPoints.col(i);
                    iter_points++;
                    if (iter_points == rnum) {
                        Neff_sampled = iter_points;
                        complete = true;
                        return;
                    }
                }
                Neff_sampled = iter_points;
                pointer = 0;
                if (rounding_requested) {
                    if (total_samples > TotalRandPoints.cols()) {
                        TotalRandPoints.conservativeResize(total_samples, P.dimension());
                        done = true;
                    }
                    TotalRandPoints.block(total_samples - min_window, 0, min_window, P.dimension()) = winPoints.transpose();
                }
                min_window = winPoints.cols();                
            }
        }
    }
};

#endif
