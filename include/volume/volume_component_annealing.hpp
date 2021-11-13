

#ifndef VOLUME_COOLING_GAUSSIANS_HPP
#define VOLUME_COOLING_GAUSSIANS_HPP


#include <iterator>
#include <vector>
#include <list>
#include <chrono>
#include <cmath>

#include "gaussian_great_cycle_walk.hpp"


template <typename NT, typename VT>
std::pair<NT, NT> get_mean_variance_vt(VT const& vec)
{
    NT mean = 0;
    NT M2 = 0;
    NT variance = 0;
    NT delta;
    NT* vec_data = vec.data();

    unsigned int M = vec.rows();
    for (int i=0; i<N; i++)
    {
        delta = (*vec_data) - mean;
        mean += delta / (i + 1);
        M2 += delta * ((*vec_data) - mean);
        variance = M2 / (i + 1);

        vec_data++;
    }
    return std::pair<NT, NT> (mean, variance);
}


template
<
    typename MT,
    typename NT,
    typename VT
>
NT get_first_two_gaussian(const MT &samples,
                          const NT &C,
                          const MT &sigma,
                          const VT &mu)
{
    NT a = NT(1), ratio;
    int N = Samples.cols();
    VT p;
    VT fn(N), fn2(N);//,NT(0.0));

    NT* fnit = fn.data();
    for (int i=0; i<N; i++)
    {
        p = samples.col(i);
        *fnit = p.dot(sigma*p));
        fnit++;
    }

    while (!done)
    {
        //a *= NT(2);
        fn2 = (-a)*fn;
        fnit = fn2.data();
        for (int i=0; i<N; i++)
        {
            *fnit = exp((*fnit));
            fnit++;
        }
        std::pair<NT, NT> mv = get_mean_variance_vt(fn2);

        // Compute a_{i+1}
        if (mv.second/(mv.first * mv.first)<=NT(1))// || mv.first/last_ratio>1.0-tol)
        {
            a *= NT(2);
        } else if (mv.second/(mv.first * mv.first)<=C)
        {
            break;
        } else
        {
            a *= NT(0.5);
        }
    }
    return a;
}



template
<
    typename GCWalk,
    typename Body,
    typename VT,
    typename MT,
    typename NT,
    typename RandomNumberGenerator
>
std::pair<NT, bool> get_next_gaussian(Body const& P,
                                      VT &p,
                                      const VT &mu,
                                      const MT &sigma,
                                      NT const& a,
                                      const unsigned int &N,
                                      const NT &ratio,
                                      const NT &C,
                                      const unsigned int& walk_length,
                                      NT &ratio_it,
                                      RandomNumberGenerator& rng)
{
    NT last_a = a;
    NT last_ratio = 0.1;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    NT k = 1.0;
    const NT tol = 0.00001;
    bool done = false;
    std::vector<NT> fn(N,NT(0.0));
    std::list<Point> randPoints;
    typedef typename std::vector<NT>::iterator viterator;
    unsigned int d = P.dimension();

    typedef typename GCWalk::template Walk
            <
                Body,
                RNGType,
                Eigen::LLT<MT>
            > CGwalk;
    
    CGwalk walk(P, p, mu, sigma, a, rng);
    MT samples(d, N);
    VT fn(N), fn2(N);

    NT* fnit = fn.data();
    for (int i=0; i<N; i++)
    {
        walk.template apply(P, p, walk_length, rng);
        samples.col(i) = p;
        *fnit = p.dot(sigma*p));
        fnit++;
    }

    bool is_last = walk.template is_inside();

    while (!done)
    {
        NT new_a = last_a * std::pow(ratio,k);
        fn2 = (-last_a/new_a) * fn;

        fnit = fn2.data();
        for (int i=0; i<N; i++)
        {
            *fnit = exp((*fnit));
            fnit++;
        }
        std::pair<NT, NT> mv = get_mean_variance_vt(fn);

        // Compute a_{i+1}
        if (mv.second/(mv.first * mv.first)>=C || mv.first/last_ratio<1.0-tol)
        {
            if (k != 1.0)
            {
                k = k / 2;
            }
            done = true;
        } else {
            k = 2 * k;
        }
        ratio_it = mv.first;
        last_ratio = mv.first;
    }
    return std::pair<NT, bool> (last_a * std::pow(ratio, k), is_last);
    //return last_a * std::pow(ratio, k);
}


template
<
    typename WalkType,
    typename Body,
    typename VT,
    typename MT,
    typename NT,
    typename RandomNumberGenerator
>
void compute_annealing_schedule(Body const& P,
                                VT &p,
                                const VT &mu,
                                const MT &sigma,
                                MT &samples,
                                NT const& ratio,
                                NT const& C,
                                NT const& frac,
                                unsigned int const& N,
                                unsigned int const& walk_length,
                                NT const& error,
                                std::vector<NT>& a_vals,
                                RandomNumberGenerator& rng)
{
    std::pair<NT, bool> res;
    std::vector<NT> ratios;
    // Compute the first gaussian
    NT a2 = get_first_two_gaussian(samples, C, sigma, mu);

    NT a1 = 0.0, a_next, ratio_it;
    const NT tol = 0.001;
    unsigned int it = 0;
    unsigned int n = P.dimension();
    //const unsigned int totalSteps = ((int)150/((1.0 - frac) * error))+1;

    //if (a_vals[0]<a_stop) a_vals[0] = a_stop;

    a_vals.push_back(a1);
    a_vals.push_back(a2);

    while (true)
    {
        // Compute the next gaussian
        res = get_next_gaussian<WalkType>(P, p, mu, sigma, a, N, ratio, C, walk_length, ratio_it, rng);
        ratios.push_back(ratio_it);
        a_vals.push_back(res.first);

        if (res.second)
        {
            break;
        }
    }
}


template <typename NT>
struct gaussian_annealing_parameters
{
    gaussian_annealing_parameters(unsigned int d)
        :   frac(0.1)
        ,   ratio(NT(1)+NT(1)/(sqrt(NT(d))))
        ,   C(NT(3))
        ,   N(500 * ((int) C) + ((int) (d * d / 2)))
        ,   W(4*d*d+500)
    {}

    NT frac;
    NT ratio;
    NT C;
    unsigned int N;
    unsigned int W;
};

template
<
    typename WalkType,
    typename MT,
    typename Body,
    typename VT,
    typename NT,
    typename RandomNumberGenerator

>
NT volume_cooling_gaussians(Body const& P,
                            VT &p,
                            VT mu,
                            MT sigma,
                            RandomNumberGenerator& rng,
                            NT const& error = 0.1,
                            unsigned int const& walk_length = 1)
{
    //const NT maxNT = std::numeric_limits<NT>::max();//1.79769e+308;
    //const NT minNT = std::numeric_limits<NT>::min();//-1.79769e+308;

    //auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    //unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    RandomNumberGenerator rng(n);

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    NT eval;
    unsigned int N = parameters.N;

    compute_annealing_schedule<WalkType>(P, p, mu, sigma, samples, ratio, C, frac, N, walk_length, error, a_vals, rng);

//#ifdef VOLESTI_DEBUG
    std::cout<<"All the variances of schedule_annealing computed in = "
            << (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    auto j=0;
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++)
    {
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;
//#endif

    // Initialization for the approximation of the ratios
    unsigned int W = parameters.W;
    unsigned int mm = a_vals.size()-1;
    std::vector<NT> last_W2(W,0);
    std::vector<NT> fn(mm,0);
    std::vector<NT> its(mm,0);
    //VT lamdas;
    //lamdas.setZero(m);
    NT vol = NT(1);
    //Point p(n); // The origin is the Chebychev center of the Polytope
    unsigned int i=0;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    viterator avalsIt = a_vals.begin();
    viterator minmaxIt;

//#ifdef VOLESTI_DEBUG
    //std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
//#endif

    //iterate over the number of ratios
    for (viterator fnIt = fn.begin();
         fnIt != fn.end();
         fnIt++, itsIt++, avalsIt++, i++)
    {
        //initialize convergence test
        bool done = false;
        NT curr_eps = error/std::sqrt((NT(mm)));
        NT min_val = std::numeric_limits<NT>::min();
        NT max_val = std::numeric_limits<NT>::max();
        unsigned int min_index = W-1;
        unsigned int max_index = W-1;
        unsigned int index = 0;
        unsigned int min_steps = 10000000;
        std::vector<NT> last_W = last_W2;

        p = mu;

        // Set the radius for the ball walk
       typedef typename GCWalk::template Walk
        <
            Body,
            RNGType,
            Eigen::LLT<MT>
        > CGwalk;
    
        CGwalk walk(P, p, mu, sigma, *avalsIt, rng);

        //update_delta<WalkType>
        //        ::apply(walk, 4.0 * radius
        //                 / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n)));

        while (!done && (*itsIt)<min_steps)
        {
            walk.template apply(P, p, *avalsIt, walk_length, rng);
            eval = p.dot(sigma*p);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + exp((-(*avalsIt)/(*(avalsIt+1)))*eval);
            NT val = (*fnIt) / (*itsIt);

            last_W[index] = val;
            if (val <= min_val)
            {
                min_val = val;
                min_index = index;
            } else if (min_index == index)
            {
                minmaxIt = std::min_element(last_W.begin(), last_W.end());
                min_val = *minmaxIt;
                min_index = std::distance(last_W.begin(), minmaxIt);
            }

            if (val >= max_val)
            {
                max_val = val;
                max_index = index;
            } else if (max_index == index)
            {
                minmaxIt = std::max_element(last_W.begin(), last_W.end());
                max_val = *minmaxIt;
                max_index = std::distance(last_W.begin(), minmaxIt);
            }

            if ( (max_val-min_val)/max_val <= curr_eps/2.0 )
            {
                done=true;
            }

            index = index%W + 1;
            if (index == W) index = 0;
        }
//#ifdef VOLESTI_DEBUG
        std::cout << "ratio " << i << " = " << (*fnIt) / (*itsIt)
                  << " N_" << i << " = " << *itsIt << std::endl;
//#endif
        vol *= ((*fnIt) / (*itsIt));
    }

//#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
//#endif

    return vol;
}



#endif
