

#ifndef VOLUME_FISCHER_ANNEALING_FAST_HPP
#define VOLUME_FISCHER_ANNEALING_FAST_HPP


#include <iterator>
#include <vector>
#include <list>
#include <chrono>
#include <cmath>

//#include "gaussian_great_cycle_walk.hpp"


template <typename NT, typename VT>
std::pair<NT, NT> get_mean_variance_fischer_vt(VT const& vec)
{
    NT mean = 0;
    NT M2 = 0;
    NT variance = 0;
    NT delta;
    const NT* vec_data = vec.data();

    unsigned int M = vec.rows();
    for (int i=0; i<M; i++)
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
NT get_first_two_fischer(const MT &samples,
                          const NT &C,
                          const NT &Cmin,
                          const VT &mu,
                          NT &ratio_it)
{
    NT a = NT(1), ratio, a1, a2;
    int N = samples.cols();
    bool done = false;
    VT p;
    VT fn(N), fn2(N);//,NT(0.0));

    NT* fnit = fn.data();
    for (int i=0; i<N; i++)
    {
        p = samples.col(i) - mu;
        *fnit = p.dot(p);
        fnit++;
    }
    int counter = 0;
    while (!done)
    {
        //counter++;
        //if (counter >100){
        //    exit(-1);
        //}
        //a *= NT(2);
        fn2 = (-a)*fn;
        fnit = fn2.data();
        for (int i=0; i<N; i++)
        {
            *fnit = exp((*fnit));
            fnit++;
        }
        //std::cout<<"fn2 = "<<fn2.transpose()<<"\n"<<std::endl;
        std::pair<NT, NT> mv = get_mean_variance_fischer_vt<NT>(fn2);
        std::cout<<"[1] a = "<<a<<std::endl;
        std::cout<<"[1] mean = "<<mv.first<<", variance = "<<mv.second<<std::endl;

        // Compute a_{i+1}
        std::cout<<"[1] var / m^2 = "<<mv.second/(mv.first * mv.first)<<std::endl;
        if ((mv.second/(mv.first * mv.first)<=C && mv.second/(mv.first * mv.first)>=Cmin) || (mv.first < 0.05))// || mv.first/last_ratio>1.0-tol)
        {
            done = true;
            ratio_it = mv.first;
            break;
        } else if (mv.second/(mv.first * mv.first)<Cmin)
        {
            a *= NT(2);
        } else
        {
            break;
        }
        
    }
    if (done)
    {
        return a;
    }
    a1 = a*0.5, a2 = a;

    while(true)
    {
        a = (a1+a2)*0.5;

        fn2 = (-a)*fn;
        fnit = fn2.data();
        for (int i=0; i<N; i++)
        {
            *fnit = exp((*fnit));
            fnit++;
        }
        //std::cout<<"fn2 = "<<fn2.transpose()<<"\n"<<std::endl;
        std::pair<NT, NT> mv = get_mean_variance_fischer_vt<NT>(fn2);
        std::cout<<"[1][BS] mean = "<<mv.first<<", variance = "<<mv.second<<std::endl;

        // Compute a_{i+1}
        std::cout<<"[1][BS] var / m^2 = "<<mv.second/(mv.first * mv.first)<<std::endl;
        if (mv.second/(mv.first * mv.first)<=C && mv.second/(mv.first * mv.first)>=Cmin)// || mv.first/last_ratio>1.0-tol)
        {
            ratio_it = mv.first;
            break;
        } else if (mv.second/(mv.first * mv.first)<Cmin)
        {
            a1 = a;
        } else
        {
            a2 = a;
        }
    }

    return a;
}



template
<
    typename GCWalk,
    typename MT,
    typename Body,
    typename VT,
    typename NT,
    typename RandomNumberGenerator
>
std::pair<NT, NT> get_next_fischer(Body const& P,
                                      VT &p,
                                      const VT &mu,
                                      NT const& a,
                                      const unsigned int &N,
                                      const NT &ratio,
                                      const NT &C,
                                      NT const& Cmin,
                                      const unsigned int& walk_length,
                                      NT &ratio_it,
                                      unsigned int const& W,
                                      RandomNumberGenerator& rng, 
                                      bool check_last = true)
{
    NT last_a = a;
    NT last_ratio = 10;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    NT k = 1.0;
    const NT tol = 0.00001, ratio_tol = 0.01;
    bool done = false;
    //std::vector<NT> fn(N,NT(0.0));
    //std::list<Point> randPoints;
    typedef typename std::vector<NT>::iterator viterator;
    unsigned int d = P.dimension();

    typedef typename GCWalk::template Walk
            <
                Body,
                RandomNumberGenerator,
                Eigen::LLT<MT>
            > CGwalk;
    
    CGwalk walk(P, p, mu, a, W, rng);
    //MT samples(d, N);
    VT fn(N), fn2(N), p_mu;

    NT* fnit = fn.data();
    for (int i=0; i<N; i++)
    {
        walk.template apply_with_check(P, p, NT(2)*a, walk_length, rng);
        p_mu = p - mu;
        //samples.col(i) = p;
        *fnit = p_mu.dot(p_mu);
        fnit++;
    }

    bool is_not_last = walk.template is_outside();
    NT ratio_outside = walk.template ratio_outside();
    NT new_a;

    std::cout<<"mu = "<<mu.transpose()<<std::endl;
    std::cout<<"N = "<<N<<std::endl;
    std::cout<<"ratio_outside = "<<ratio_outside<<std::endl;

    if (check_last) {
        if (ratio_outside < ratio_tol){
            return std::pair<NT, NT> (a, ratio_outside);
        }
    }
    int counter = 0;
    while (!done)
    {
        //counter++;
        //if (counter >10){
        //    exit(-1);
        //}
        new_a = last_a * std::pow(ratio,k);
        std::cout<<"new_a = "<<new_a<<std::endl;
        fn2 = (last_a - new_a) * fn;

        fnit = fn2.data();
        for (int i=0; i<N; i++)
        {
            *fnit = exp((*fnit));
            fnit++;
        }
        //std::cout<<"fn2 = "<<fn2.transpose()<<"\n"<<std::endl;
        std::pair<NT, NT> mv = get_mean_variance_fischer_vt<NT>(fn2);
        std::cout<<"mean = "<<mv.first<<", variance = "<<mv.second<<std::endl;

        // Compute a_{i+1}
        std::cout<<"var / m^2 = "<<mv.second/(mv.first * mv.first)<<std::endl;
        std::cout<<"mv.first/last_ratio = "<<mv.first/last_ratio<<std::endl;
        //exit(-1);
        //std::cout<<"C = "<<C<<std::endl;
        //std::cout<<"ratio = "<<ratio<<std::endl;

        if (mv.second/(mv.first * mv.first)>=C)// || mv.first/last_ratio>1.0-tol)
        {
            //if (k != 1.0)
            //{
            //    k = k / 2;
            //}
            ratio_it = mv.first;
            done = true;
        } else {
            k = 2 * k;
        }
        last_ratio = mv.first;
    }

    NT k1 = k/2, k2 = k;

    counter = 0;
    while(true) 
    {
        counter++;
        
        k = (k1+k2)/NT(2);
        new_a = last_a * std::pow(ratio,k);
        std::cout<<"new_a = "<<new_a<<std::endl;
        fn2 = (last_a - new_a) * fn;

        fnit = fn2.data();
        for (int i=0; i<N; i++)
        {
            *fnit = exp((*fnit));
            fnit++;
        }
        //std::cout<<"fn2 = "<<fn2.transpose()<<"\n"<<std::endl;
        std::pair<NT, NT> mv = get_mean_variance_fischer_vt<NT>(fn2);
        std::cout<<"[BS] mean = "<<mv.first<<", variance = "<<mv.second<<std::endl;

        // Compute a_{i+1}
        std::cout<<"[BS] var / m^2 = "<<mv.second/(mv.first * mv.first)<<std::endl;
        std::cout<<"[BS] mv.first/last_ratio = "<<mv.first/last_ratio<<std::endl;
        //exit(-1);
        //std::cout<<"C = "<<C<<std::endl;
        //std::cout<<"ratio = "<<ratio<<std::endl;

        if (mv.second/(mv.first * mv.first)>=Cmin && mv.second/(mv.first * mv.first)<=C)// || mv.first/last_ratio>1.0-tol)
        {
            //if (k != 1.0)
            //{
            //    k = k / 2;
            //}
            ratio_it = mv.first;
            break;
        } else if (mv.second/(mv.first * mv.first)>C) {
            k2 = k;
        } else {
            k1 = k;
        }
        if (counter >10){
            ratio_it = mv.first;
            break;
        }
        last_ratio = mv.first;
    }
    

    return std::pair<NT, NT> (last_a * std::pow(ratio, k), ratio_outside);
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
void compute_annealing_schedule_fischer(Body const& P,
                                VT &p,
                                const VT &mu,
                                MT const& samples,
                                NT const& ratio,
                                NT const& C,
                                NT const& Cmin,
                                unsigned int const& N,
                                unsigned int const& walk_length,
                                std::vector<NT>& a_vals,
                                std::vector<NT> &ratios,
                                unsigned int const& W,
                                RandomNumberGenerator& rng)
{
    std::pair<NT, bool> res;
    NT ratio_it;
    // Compute the first gaussian
    NT a2 = get_first_two_fischer(samples, C, Cmin, mu, ratio_it);
    ratios.push_back(ratio_it);
    std::cout<<"ratio = "<<ratio_it<<std::endl;
    std::cout<<"first two computed"<<"\n"<<std::endl;
    NT a1 = 0.0, a_next;
    //const NT tol = 0.001;
    unsigned int it = 1;
    unsigned int n = P.dimension();
    //const unsigned int totalSteps = ((int)150/((1.0 - frac) * error))+1;

    //if (a_vals[0]<a_stop) a_vals[0] = a_stop;

    a_vals.push_back(a1);
    a_vals.push_back(a2);
    VT p0=p;

    while (true)
    {
        // Compute the next gaussian
        p=p0;
        res = get_next_fischer<WalkType, MT>(P, p, mu, a_vals[it], N, ratio, C, Cmin, walk_length, ratio_it, W, rng);

        if (res.second < 0.01)
        {
            break;
        }
        
        ratios.push_back(ratio_it);
        a_vals.push_back(res.first);
        std::cout<<"a_next = "<<res.first<<std::endl;
        std::cout<<"ratio = "<<ratio_it<<std::endl;
        std::cout<<"num of phases = "<<a_vals.size()<<"\n"<<std::endl;
        //if(a_vals.size()>30){
         //   exit(-1);
        //}
        it++;        
    }
}


template <typename NT>
struct gaussian_annealing_parameters
{
    gaussian_annealing_parameters(unsigned int d)
        :   frac(0.1)
        ,   ratio(NT(1)+(NT(1)/NT(d)))
        ,   C(NT(2))
        ,   Cmin(NT(1))
        ,   N(500 * ((int) C) + ((int) (d * d / 2)))
        ,   W(4*d*d+1000)
    {}

    NT frac;
    NT ratio;
    NT C;
    NT Cmin;
    unsigned int N;
    unsigned int W;
};


template <typename MT, typename VT>
MT estimate_cov(std::vector<VT> &points_temp, unsigned int const& d)
{
    VT avg, temp(d);
    MT sigma;

    typedef typename std::vector<VT>::iterator viterator;
    sigma.setZero(d,d);
    avg.setZero(d);
    unsigned int N = 0;

    for (viterator fnIt = points_temp.begin();
         fnIt != points_temp.end();
         fnIt++)
    {
        avg += (*fnIt);
        N++;
    }
    avg *= (1.0/double(N));

    for (viterator fnIt = points_temp.begin();
         fnIt != points_temp.end();
         fnIt++)
    {
        temp = (*fnIt) - avg;
        sigma += temp * temp.transpose();
    }
    sigma *= (1.0/double(N));

    return sigma;
}


template
<
    typename UniformWalkType,
    typename WalkType,
    typename MT,
    typename Body,
    typename VT,
    typename NT,
    typename RandomNumberGenerator
>
std::pair<NT, NT> volume_component_cooling_fischer(Body const& P,
                                      VT &p,
                                      VT const& mu,
                                      RandomNumberGenerator& rng,
                                      unsigned int const& WW,
                                      NT const& error = 0.1,
                                      unsigned int const& walk_length = 1)
{
    //const NT maxNT = std::numeric_limits<NT>::max();//1.79769e+308;
    //const NT minNT = std::numeric_limits<NT>::min();//-1.79769e+308;

    //auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    //unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    //RandomNumberGenerator rng(n);

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    NT Cmin = parameters.Cmin;
    NT eval;
    VT p_mu;
    unsigned int N = parameters.N;
    std::vector<NT> ratios;
    MT samples(n, N);
    VT p0 = p;

    typedef typename UniformWalkType::template Walk
        <
            Body,
            RandomNumberGenerator
        > CGWalk;
    CGWalk Uwalk(P, p, rng);


    for (int jj=0; jj<N; jj++)
    {
        Uwalk.template apply(P, p, walk_length, rng);
        samples.col(jj) = p;
    }
    p = p0;

    compute_annealing_schedule_fischer<WalkType>(P, p, mu, samples, ratio, C, Cmin, N, walk_length, a_vals, ratios, WW, rng);

    int j=0, MM = samples.cols();
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++)
    {
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;

    //std::sort(a_vals.begin(), a_vals.end(), std::greater<NT>());

//#ifdef VOLESTI_DEBUG
    //std::cout<<"All the variances of schedule_annealing computed in = "
            //<< (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    j=0;
    for (auto avalIt = ratios.begin(); avalIt!=ratios.end(); avalIt++, j++)
    {
        std::cout<<"r_"<<j<<" = "<<*avalIt<<" ";
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
    unsigned int i=0, n_max = 50*n*n;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    viterator avalsIt = a_vals.begin();
    viterator minmaxIt;

//#ifdef VOLESTI_DEBUG
    //std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
//#endif

    typedef typename WalkType::template Walk
        <
            Body,
            RandomNumberGenerator,
            Eigen::LLT<MT>
        > CGwalk;

    CGwalk walk(P, p, mu, NT(2)*(*avalsIt), WW, rng);
    std::vector<VT> points;
    MT sigma_temp = MT::Identity(n,n);

    p0 = mu;
    Uwalk.template initialize(P, p0, rng);

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
        std::vector<VT> points_temp = points;
        std::vector<NT> last_W = last_W2;

        p = mu;

        // Set the radius for the ball walk
       
        walk.template set_sigma(sigma_temp);
        walk.template initialize(P, p, NT(2)*(*avalsIt), rng);

        while (!done && (*itsIt)<min_steps)
        {
            if (i==0) {
                Uwalk.template apply(P, p0, walk_length, rng);
                if ((*itsIt) <= NT(n_max)) {
                    points_temp.push_back(p0);
                }
                p_mu.noalias() = p0 - mu;
            } else {
                walk.template apply(P, p, NT(2)*(*avalsIt), walk_length, rng);
                if ((*itsIt) <= NT(n_max)) {
                    points_temp.push_back(p);
                }
                p_mu.noalias() = p - mu;
            }
            eval = p_mu.dot(p_mu);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + exp( ( (*avalsIt) - (*(avalsIt+1)) ) * eval);
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
                  << " N_" << i << " = " << *itsIt << ", mm = " << mm << std::endl;
//#endif
        vol *= ((*fnIt) / (*itsIt));
        if (i<mm-1) {
            sigma_temp = estimate_cov<MT>(points_temp, n);
        }
    }

//#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
//#endif

    return std::pair<NT, NT>(vol, a_vals[a_vals.size()-1]);
}



template
<
    typename UniformWalkType,
    typename WalkType,
    typename MT,
    typename NT,
    typename Body,
    typename VT,
    typename RandomNumberGenerator
>
std::pair<VT, VT> compute_annealing_fischer(Body const& P,
                                      VT &p,
                                      VT const& mu,
                                      RandomNumberGenerator& rng,
                                      unsigned int const& WW,
                                      unsigned int const& walk_length = 1)
{
    //const NT maxNT = std::numeric_limits<NT>::max();//1.79769e+308;
    //const NT minNT = std::numeric_limits<NT>::min();//-1.79769e+308;

    //auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    //unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    //RandomNumberGenerator rng(n);

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    NT Cmin = parameters.Cmin;
    unsigned int N = parameters.N;
    std::vector<NT> ratios;
    MT samples(n, N);
    VT p0=p;

    typedef typename UniformWalkType::template Walk
        <
            Body,
            RandomNumberGenerator
        > CGWalk;
    CGWalk walk(P, p, rng);

    for (int jj=0; jj<N; jj++)
    {
        walk.template apply(P, p, walk_length, rng);
        samples.col(jj) = p;
    }
    p = p0;
    compute_annealing_schedule_fischer<WalkType>(P, p, mu, samples, ratio, C, Cmin, N, walk_length, a_vals, ratios, WW, rng);

    VT a_sequence(a_vals.size());
    int j=0;
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++)
    {
        a_sequence(j) = (*avalIt);
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;

    //std::sort(a_vals.begin(), a_vals.end(), std::greater<NT>());

//#ifdef VOLESTI_DEBUG
    //std::cout<<"All the variances of schedule_annealing computed in = "
            //<< (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    VT ratios_vt(ratios.size());
    j=0;
    for (auto avalIt = ratios.begin(); avalIt!=ratios.end(); avalIt++, j++)
    {
        ratios_vt(j) = (*avalIt);
        std::cout<<"r_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;
//#endif

    return std::pair<VT, VT>(a_sequence, ratios_vt);
}


template
<
    typename UniformWalkType,
    typename WalkType,
    typename MT,
    typename Body,
    typename VT,
    typename NT,
    typename RandomNumberGenerator
>
NT estimate_ratios_fischer(Body const& P,
                                          VT &p,
                                          VT const& mu,
                                          VT const& a_vals,
                                          VT &ratios,
                                          unsigned int& N,
                                          RandomNumberGenerator& rng,
                                          unsigned int const& WW,
                                          NT const& error = 0.1,
                                          unsigned int const& walk_length = 1)
{
    unsigned int n = P.dimension();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    // Initialization for the approximation of the ratios
    unsigned int W = parameters.W;
    unsigned int mm = a_vals.rows()-1;
    std::vector<NT> last_W2(W,0);
    
    std::vector<NT> its(mm,0);
    ratios *= NT(N);
    VT p_mu;
    std::vector<NT> fn(mm,0);
    //VT lamdas;
    //lamdas.setZero(m);
    NT vol = NT(1), eval;
    //Point p(n); // The origin is the Chebychev center of the Polytope
    unsigned int i=0, count, n_max = 50*n*n;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    const NT* avalsIt = a_vals.data();
    //NT* fnIt = fn.data();
    viterator minmaxIt;

//#ifdef VOLESTI_DEBUG
    //std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
//#endif

    typedef typename WalkType::template Walk
        <
            Body,
            RandomNumberGenerator,
            Eigen::LLT<MT>
        > CGwalk;

    CGwalk walk(P, p, mu, NT(2)*(*avalsIt), WW, rng);

    VT p0 = mu;
    typedef typename UniformWalkType::template Walk
        <
            Body,
            RandomNumberGenerator
        > CGWalk;
    CGWalk Uwalk(P, p0, rng);

    std::vector<VT> points;
    MT sigma_temp = MT::Identity(n, n);

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
        std::vector<VT> points_temp = points;
        std::vector<NT> last_W = last_W2;

        p = mu;

        // Set the radius for the ball walk
        
        walk.template set_sigma(sigma_temp);
        walk.template initialize(P, p, NT(2)*(*avalsIt), rng);

        while (!done && (*itsIt)<min_steps)
        {
            if (i==0) {
                Uwalk.template apply(P, p0, walk_length, rng);
                if ((*itsIt) <= NT(n_max)) {
                    points_temp.push_back(p0);
                }
                p_mu.noalias() = p0 - mu;
            } else {
                walk.template apply(P, p, NT(2)*(*avalsIt), walk_length, rng);
                if ((*itsIt) <= NT(n_max)) {
                    points_temp.push_back(p);
                }
                p_mu.noalias() = p - mu;
            }
            eval = p_mu.dot(p_mu);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + exp( ( (*avalsIt) - (*(avalsIt+1)) ) * eval);
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
                  << " N_" << i << " = " << *itsIt << ", mm = " << mm << std::endl;
//#endif
        vol *= (*fnIt) / (*itsIt);
        //vol *= ((*fnIt) / (*itsIt));
        if (i<mm-1) {
            sigma_temp = estimate_cov<MT>(points_temp, n);
        }
    }

//#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
        //std::cout<<"\nvol = "<<vol<<"\n"<<std::endl;
//#endif

    return vol;
}






template
<
    typename WalkType,
    typename MT,
    typename Body,
    typename VT,
    typename NT,
    typename RandomNumberGenerator
>
std::pair<NT, NT> related_volume_cooling_fischer(Body const& P,
                                      VT &p,
                                      VT const& mu,
                                      NT const& a_max,
                                      NT const& a_min,
                                      NT &ratio_min,
                                      NT &ratio_max,
                                      RandomNumberGenerator& rng,
                                      unsigned int const& WW,
                                      NT const& error = 0.1,
                                      unsigned int const& walk_length = 1)
{
    //const NT maxNT = std::numeric_limits<NT>::max();//1.79769e+308;
    //const NT minNT = std::numeric_limits<NT>::min();//-1.79769e+308;

    //auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension(), it = 0;
    //unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());
    //RandomNumberGenerator rng(n);

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    NT Cmin = parameters.Cmin;
    NT eval, ratio_it;
    VT p_mu;
    unsigned int N = parameters.N;
    std::vector<NT> ratios;

    //compute_annealing_schedule_fischer<WalkType>(P, p, mu, samples, ratio, C, Cmin, N, walk_length, a_vals, ratios, WW, rng);

    a_vals.push_back(a_min);
    //a_vals.push_back(a2);
    VT p0 = p;
    NT ratio_min_temp = NT(1);
    NT ratio_max_temp = ratio_max;
    std::pair<NT,bool> res;

    while (true)
    {
        // Compute the next gaussian
        p = p0;
        res = get_next_fischer<WalkType, MT>(P, p, mu, a_vals[it], N, ratio, C, Cmin, walk_length, ratio_it, WW, rng, false);
        
        ratios.push_back(ratio_it);
        ratio_min_temp *= (NT(1) / ratio_it);

        std::cout<<"a_next = "<<res.first<<std::endl;
        std::cout<<"ratio = "<<ratio_it<<std::endl;
        std::cout<<"num of phases = "<<a_vals.size()<<"\n"<<std::endl;

        if (res.first < a_max) {
            a_vals.push_back(res.first);
        } else {
            a_vals.push_back(a_max);
            break;
        }

        it++;        
    }

    int j=0;
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++)
    {
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;

    j=0;
    for (auto avalIt = ratios.begin(); avalIt!=ratios.end(); avalIt++, j++)
    {
        std::cout<<"r_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;

    if ((ratio_min_temp*ratio_min) / ratio_max_temp < 1e-05 || (ratio_min_temp*ratio_min) / ratio_max_temp > 1e05){
        std::cout<<"ratio_min_temp / ratio_max_temp = "<<(ratio_min_temp*ratio_min) / ratio_max_temp<<std::endl;
        return std::pair<NT, NT> (ratio_min_temp, ratio_max_temp);
    } 
    ratio_max_temp = ratio_max;
    ratio_min_temp = NT(1);

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

    typedef typename WalkType::template Walk
        <
            Body,
            RandomNumberGenerator,
            Eigen::LLT<MT>
        > CGwalk;

    CGwalk walk(P, p, mu, NT(2)*(*avalsIt), WW, rng);
    std::vector<VT> points;

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
        //std::vector<VT> points_temp = points;
        std::vector<NT> last_W = last_W2;

        p = mu;

        // Set the radius for the ball walk
       
    
        //walk.template set_sigma(sigma_temp);
        walk.template initialize(P, p, NT(2)*(*avalsIt), rng);

        //update_delta<WalkType>
        //        ::apply(walk, 4.0 * radius
        //                 / std::sqrt(std::max(NT(1.0), *avalsIt) * NT(n)));

        while (!done && (*itsIt)<min_steps)
        {
            walk.template apply_ratio_esti(P, p, NT(2)*(*avalsIt), walk_length, rng);
            //points_temp.push_back(p);
            p_mu.noalias() = p - mu;
            eval = p_mu.dot(p_mu);

            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + exp( ( (*avalsIt) - (*(avalsIt+1)) ) * eval);
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
                  << " N_" << i << " = " << *itsIt << ", mm = " << mm << std::endl;
//#endif
        //vol *= ((*fnIt) / (*itsIt));
        ratio_min_temp *= ((*fnIt) / (*itsIt));
        //if (i<mm-1) {
        //    sigma_temp = estimate_cov<MT>(points_temp, n);
        //}
    }

//#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
//#endif

    return std::pair<NT, NT> (ratio_min_temp, ratio_max_temp);
}

#endif
