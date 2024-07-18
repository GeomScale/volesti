// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
// Copyright (c) 2024 Vladimir Necula

// Contributed and/or modified by Vladimir Necula, as part of Google Summer of
// Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with
// Riemannian Hamiltonian
// Monte Carlo in a Constrained Space"
#ifndef VOLUME_COOLING_GAUSSIANS_CRHMC_HPP
#define VOLUME_COOLING_GAUSSIANS_CRHMC_HPP

//#define VOLESTI_DEBUG

#include "volume/volume_cooling_gaussians.hpp"
#include "preprocess/crhmc/crhmc_problem.h"

template
<
    typename CRHMCWalkType,
    typename crhmc_walk_params,
    int simdLen,
    typename Grad,
    typename Func,
    typename CrhmcProblem,
    typename Polytope,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
NT get_next_gaussian(Polytope& P,
                     Point &p,
                     NT const& a,
                     const unsigned int &N,
                     const NT &ratio,
                     const NT &C,
                     const unsigned int& walk_length,
                     RandomNumberGenerator& rng,
                     Grad& g,
                     Func& f,
                     crhmc_walk_params& parameters,
                     CrhmcProblem& problem,
                     CRHMCWalkType& crhmc_walk)
{

    NT last_a = a;
    NT last_ratio = 0.1;
    //k is needed for the computation of the next variance a_{i+1} = a_i * (1-1/d)^k
    NT k = 1.0;
    const NT tol = 0.00001;
    bool done=false;
    std::vector<NT> fn(N,NT(0.0));
    std::list<Point> randPoints;
    typedef typename std::vector<NT>::iterator viterator;

    //sample N points
    PushBackWalkPolicy push_back_policy;
    bool raw_output = false; 

    typedef CrhmcRandomPointGenerator<CRHMCWalkType> CRHMCRandomPointGenerator;

    CRHMCRandomPointGenerator::apply(problem, p, N, walk_length, randPoints,
                                    push_back_policy, rng, g, f, parameters, crhmc_walk, simdLen, raw_output);

    while (!done)
    {
        NT new_a = last_a * std::pow(ratio,k);

        auto fnit = fn.begin();
        for (auto pit=randPoints.begin(); pit!=randPoints.end(); ++pit, fnit++)
        {
            *fnit = eval_exp(*pit, new_a)/eval_exp(*pit, last_a);
        }
        std::pair<NT, NT> mv = get_mean_variance(fn);

        // Compute a_{i+1}
        if (mv.second/(mv.first * mv.first)>=C || mv.first/last_ratio<1.0+tol)
        {
            if (k != 1.0)
            {
                k = k / 2;
            }
            done = true;
        } 
        else {
            k = 2 * k;
        }
        last_ratio = mv.first;
    }
    return last_a * std::pow(ratio, k);
}

// Compute the sequence of spherical gaussians
template
<
    int simdLen,
    typename Polytope,
    typename NT,
    typename RandomNumberGenerator
>
void compute_annealing_schedule(Polytope& P,
                                NT const& ratio,
                                NT const& C,
                                NT const& frac,
                                unsigned int const& N,
                                unsigned int const& walk_length,
                                NT const& chebychev_radius,
                                NT const& error,
                                std::vector<NT>& a_vals,
                                RandomNumberGenerator& rng)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::MT 	MT;
    
    typedef typename GaussianFunctor::FunctionFunctor<Point>    Func;
    typedef typename GaussianFunctor::GradientFunctor<Point>    Grad;
    typedef typename GaussianFunctor::HessianFunctor<Point>     Hess;
    typedef typename GaussianFunctor::parameters<NT, Point>     func_params;

    typedef crhmc_input<MT, Point, Func, Grad, Hess> Input;
    typedef crhmc_problem<Point, Input> CrhmcProblem;   

    typedef ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad, simdLen> Solver;

    typedef typename CRHMCWalk::template Walk
            <
                    Point,
                    CrhmcProblem,
                    RandomNumberGenerator,
                    Grad,
                    Func,
                    Solver
            > CRHMCWalkType;

    typedef typename CRHMCWalk::template parameters
          <
                  NT,
                  Grad
          > crhmc_walk_params;

    typedef CrhmcRandomPointGenerator<CRHMCWalkType> CRHMCRandomPointGenerator;


    // Compute the first gaussian
    // This uses the function from the standard volume_cooling_gaussians.hpp
    get_first_gaussian(P, frac, chebychev_radius, error, a_vals);
    NT a_stop = 0.0;
    const NT tol = 0.001;
    unsigned int it = 0;
    unsigned int n = P.dimension();
    const unsigned int totalSteps = ((int)150/((1.0 - frac) * error))+1;

    if (a_vals[0]<a_stop) a_vals[0] = a_stop;

#ifdef VOLESTI_DEBUG
    std::cout<<"first gaussian computed "<< a_vals[0] <<std::endl;
#endif

#ifdef VOLESTI_DEBUG
    std::cout<<"Computing the sequence of gaussians..\n"<<std::endl;
#endif

    while (true)
    {

        NT curr_fn = 0;
        NT curr_its = 0;
        auto steps = totalSteps;

        //TODO: potential problem creation and preprocessing optimization
        
        //create the crhmc problem for this variance
        int dimension = P.dimension();
        func_params f_params = func_params(Point(dimension), a_vals[it], 1);
        
        Func f(f_params);
        Grad g(f_params);
        Hess h(f_params);

        Input input = convert2crhmc_input<Input, Polytope, Func, Grad, Hess>(P, f, g, h);

        typedef crhmc_problem<Point, Input> CrhmcProblem;
        CrhmcProblem problem = CrhmcProblem(input);

        Point p = Point(problem.center);

        if(problem.terminate){return;}

        problem.options.simdLen = simdLen;
        crhmc_walk_params params(input.df, p.dimension(), problem.options);

        if (input.df.params.eta > 0) {
            params.eta = input.df.params.eta;
        }

        int dim;
        dim = p.dimension();

        //create the walk object for this problem
        CRHMCWalkType walk = CRHMCWalkType(problem, p, input.df, input.f, params);

        // Compute the next gaussian
        NT next_a = get_next_gaussian<CRHMCWalkType, crhmc_walk_params, simdLen, Grad, Func, CrhmcProblem>
                      (P, p, a_vals[it], N, ratio, C, walk_length, rng, g, f, params, problem, walk);

#ifdef VOLESTI_DEBUG
    std::cout<<"Next Gaussian " << next_a <<std::endl;
#endif

        // Compute some ratios to decide if this is the last gaussian
        for (unsigned  int j = 0; j < steps; j++)
        {
            walk.apply(rng, walk_length);
            p = walk.getPoint();
#ifdef VOLESTI_DEBUG
    std::cout<<"Walk point " << std::endl;
    p.print();
#endif
            curr_its += 1.0;
            curr_fn += eval_exp(p, next_a) / eval_exp(p, a_vals[it]);
        }

#ifdef VOLESTI_DEBUG
    std::cout<<"Condition function ratio " << curr_fn/curr_its <<std::endl;
#endif

        // Remove the last gaussian.
        // Set the last a_i equal to zero
        if (next_a>0 && curr_fn/curr_its>(1.0+tol))
        {
            a_vals.push_back(next_a);
            it++;
        } else if (next_a <= 0)
        {
            a_vals.push_back(a_stop);
            it++;
            break;
        } else {
            a_vals[it] = a_stop;
            break;
        }
    }
#ifdef VOLESTI_DEBUG
    std::cout<<"first gaussian after WHILE "<< a_vals[0] <<std::endl;
#endif
}

template
<
    typename Polytope,
    typename RandomNumberGenerator,
    int simdLen = 8
>
double volume_cooling_gaussians(Polytope& Pin,
                                RandomNumberGenerator& rng,
                                double const& error = 0.1,
                                unsigned int const& walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT 	NT;
    typedef typename Polytope::VT 	VT;
    typedef typename Polytope::MT 	MT;
    typedef typename GaussianFunctor::FunctionFunctor<Point>    Func;
    typedef typename GaussianFunctor::GradientFunctor<Point>    Grad;
    typedef typename GaussianFunctor::HessianFunctor<Point>     Hess;
    typedef typename GaussianFunctor::parameters<NT, Point>     func_params;

    typedef crhmc_input<MT, Point, Func, Grad, Hess> Input;
    typedef crhmc_problem<Point, Input> CrhmcProblem;   

    typedef ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad, simdLen> Solver;

    typedef typename CRHMCWalk::template Walk
            <
                    Point,
                    CrhmcProblem,
                    RandomNumberGenerator,
                    Grad,
                    Func,
                    Solver
            > CRHMCWalkType;

    typedef typename CRHMCWalk::template parameters
          <
                  NT,
                  Grad
          > crhmc_walk_params;
    
    typedef CrhmcRandomPointGenerator<CRHMCWalkType> CRHMCRandomPointGenerator;

    auto P(Pin); //copy and work with P because we are going to shift
    unsigned int n = P.dimension();
    unsigned int m = P.num_of_hyperplanes();
    gaussian_annealing_parameters<NT> parameters(P.dimension());

    // Consider Chebychev center as an internal point
    auto InnerBall = P.ComputeInnerBall();
    if (InnerBall.second < 0.0) return -1.0;

    Point c = InnerBall.first;
    NT radius = InnerBall.second;

    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(c.getCoefficients());

    // Computing the sequence of gaussians
#ifdef VOLESTI_DEBUG
    std::cout<<"\n\nComputing annealing...\n"<<std::endl;
    double tstart2 = (double)clock()/(double)CLOCKS_PER_SEC;
#endif

    // Initialization for the schedule annealing
    std::vector<NT> a_vals;
    NT ratio = parameters.ratio;
    NT C = parameters.C;
    unsigned int N = parameters.N;

    compute_annealing_schedule<simdLen>(P, ratio, C, parameters.frac, N, walk_length, radius, error, a_vals, rng);

#ifdef VOLESTI_DEBUG
    std::cout<<"All the variances of schedule_annealing computed in = "
            << (double)clock()/(double)CLOCKS_PER_SEC-tstart2<<" sec"<<std::endl;
    auto j=0;
    for (auto avalIt = a_vals.begin(); avalIt!=a_vals.end(); avalIt++, j++){
        std::cout<<"a_"<<j<<" = "<<*avalIt<<" ";
    }
    std::cout<<std::endl<<std::endl;
#endif

    // Initialization for the approximation of the ratios
    unsigned int W = parameters.W;
    unsigned int mm = a_vals.size()-1;
    std::vector<NT> last_W2(W,0);
    std::vector<NT> fn(mm,0);
    std::vector<NT> its(mm,0);
    VT lamdas;
    lamdas.setZero(m);
    NT vol = std::pow(M_PI/a_vals[0], (NT(n))/2.0);
    unsigned int i=0;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    viterator avalsIt = a_vals.begin();
    viterator minmaxIt;


#ifdef VOLESTI_DEBUG
    std::cout<<"volume of the first gaussian = "<<vol<<"\n"<<std::endl;
    std::cout<<"computing ratios..\n"<<std::endl;
#endif

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
        unsigned int min_steps = 0;
        std::vector<NT> last_W = last_W2;

        // Set the radius for the ball walk
        //creating the walk object
        int dimension = P.dimension();
        func_params f_params = func_params(Point(dimension), *avalsIt, 1);
        
        Func f(f_params);
        Grad g(f_params);
        Hess h(f_params);

        //create the crhmc problem
        Input input = convert2crhmc_input<Input, Polytope, Func, Grad, Hess>(P, f, g, h);

        CrhmcProblem problem = CrhmcProblem(input);

        Point p = Point(problem.center);

        if(problem.terminate){return 0;}

        problem.options.simdLen=simdLen;
        crhmc_walk_params params(input.df, p.dimension(), problem.options);

        if (input.df.params.eta > 0) {
            params.eta = input.df.params.eta;
        }

        CRHMCWalkType walk = CRHMCWalkType(problem, p, input.df, input.f, params);

        while (!done || (*itsIt)<min_steps)
        {
            walk.apply(rng, walk_length);
            p = walk.getPoint();
            *itsIt = *itsIt + 1.0;
            *fnIt = *fnIt + eval_exp(p,*(avalsIt+1)) / eval_exp(p,*avalsIt);
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
#ifdef VOLESTI_DEBUG
        std::cout << "ratio " << i << " = " << (*fnIt) / (*itsIt)
                  << " N_" << i << " = " << *itsIt << std::endl;
#endif
        vol *= ((*fnIt) / (*itsIt));
    }

#ifdef VOLESTI_DEBUG
        NT sum_of_steps = 0.0;
        for(viterator it = its.begin(); it != its.end(); ++it) {
            sum_of_steps += *it;
        }
        auto steps= int(sum_of_steps);
        std::cout<<"\nTotal number of steps = "<<steps<<"\n"<<std::endl;
#endif

    return vol;
}

#endif // VOLUME_COOLING_GAUSSIANS_CRHMC_HPP
