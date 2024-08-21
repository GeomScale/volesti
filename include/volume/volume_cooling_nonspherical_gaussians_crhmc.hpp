// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2024 Vissarion Fisikopoulos
// Copyright (c) 2018-2024 Apostolos Chalkis
// Copyright (c) 2024 Vladimir Necula

// Contributed and/or modified by Vladimir Necula, as part of Google Summer of
// Code 2024 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLUME_COOLING_NON_SPHERICAL_GAUSSIANS_CRHMC_HPP
#define VOLUME_COOLING_NON_SPHERICAL_GAUSSIANS_CRHMC_HPP

//#define VOLESTI_DEBUG

#include "volume/volume_cooling_gaussians.hpp"
#include "preprocess/crhmc/crhmc_problem.h"
#include "preprocess/max_inscribed_ellipsoid.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"

////////////////////////////// Algorithms

// Gaussian Anealling

// Compute a_{i+1} when a_i is given

template
<
    typename CRHMCWalkType,
    typename crhmc_walk_params,
    int simdLen,
    typename Grad,
    typename Func,
    typename CrhmcProblem,
    typename Point,
    typename NT,
    typename RandomNumberGenerator
>
void burnIn(Point &p,
            const unsigned int& walk_length,
            RandomNumberGenerator& rng,
            Grad& g,
            Func& f,
            crhmc_walk_params& parameters,
            CrhmcProblem& problem,
            CRHMCWalkType& crhmc_walk,
            NT burnIn_sample) 
{
    std::list<Point> randPoints;

    // burnIn
    PushBackWalkPolicy push_back_policy;
    bool raw_output = false;
    typedef CrhmcRandomPointGenerator<CRHMCWalkType> CRHMCRandomPointGenerator;

    CRHMCRandomPointGenerator::apply(problem, p, burnIn_sample, walk_length, randPoints,
                                    push_back_policy, rng, g, f, parameters, crhmc_walk, simdLen, raw_output);

}


template
<
    int simdLen,
    typename Polytope,
    typename Point,
    typename NT,
    typename MT,
    typename RandomNumberGenerator
>
NT get_next_gaussian(Polytope& P,
                    Point& p,
                    NT const& a,
                    const unsigned int &N,
                    const NT &ratio,
                    const NT &C,
                    const unsigned int& walk_length,
                    RandomNumberGenerator& rng,
                    MT const& inv_covariance_matrix) 
{
    typedef typename NonSphericalGaussianFunctor::FunctionFunctor<Point>    Func;
    typedef typename NonSphericalGaussianFunctor::GradientFunctor<Point>    Grad;
    typedef typename NonSphericalGaussianFunctor::HessianFunctor<Point>     Hess;
    typedef typename NonSphericalGaussianFunctor::parameters<NT, Point>     func_params;

    typedef crhmc_input<MT, Point, Func, Grad, ZeroFunctor<Point>> Input;

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

    // Create the CRHMC problem for this variance
    int dimension = P.dimension();
    func_params f_params = func_params(Point(dimension), a, -1, inv_covariance_matrix);

    Func f(f_params);
    Grad g(f_params);
    Hess h(f_params);
    ZeroFunctor<Point> zerof;

    Input input = convert2crhmc_input<Input, Polytope, Func, Grad, ZeroFunctor<Point>>(P, f, g, zerof);

    typedef crhmc_problem<Point, Input> CrhmcProblem;
    
    CrhmcProblem problem = CrhmcProblem(input);

    if(problem.terminate) { return 0;}

    problem.options.simdLen = simdLen;

    crhmc_walk_params params(input.df, p.dimension(), problem.options);

    // Create the walk object for this problem
    CRHMCWalkType walk = CRHMCWalkType(problem, p, input.df, input.f, params);

    burnIn<CRHMCWalkType, crhmc_walk_params, simdLen, Grad, Func, CrhmcProblem>(
        p, walk_length, rng, g, f, params, problem, walk, 2.5 * N);

    NT last_a = a;
    NT last_ratio = 0.1;
    NT k = 1.0;
    const NT tol = 0.00001;
    bool done = false;
    std::vector<NT> fn(N, NT(0.0));
    std::list<Point> randPoints;

    // sample N points
    PushBackWalkPolicy push_back_policy;
    bool raw_output = false;
    typedef CrhmcRandomPointGenerator<CRHMCWalkType> CRHMCRandomPointGenerator;

    CRHMCRandomPointGenerator::apply(problem, p, N, walk_length, randPoints,
                                    push_back_policy, rng, g, f, params, walk, simdLen, raw_output);

    while (!done) {
        NT new_a = last_a * std::pow(ratio, k);

        auto fnit = fn.begin();
        for (auto pit = randPoints.begin(); pit != randPoints.end(); ++pit, fnit++) {
            *fnit = eval_exp(*pit, inv_covariance_matrix, new_a, last_a);
        }

        std::pair<NT, NT> mv = get_mean_variance(fn);

        // Compute a_{i+1}
        if (mv.second / (mv.first * mv.first) >= C || mv.first / last_ratio < 1.0 + tol) {
            if (k != 1.0) {
                k = k / 2;
            }
            done = true;
        } else {
            k = 2 * k;
        }

        last_ratio = mv.first;
    }

    // Return the new a value as a scalar
    return last_a * std::pow(ratio, k);
}

// Compute the sequence of non spherical gaussians
template<
    int simdLen,
    typename Point,
    typename Polytope,
    typename NT,
    typename MT,
    typename RandomNumberGenerator
>
void compute_annealing_schedule(Polytope Pin_copy,
                                Polytope P_copy,
                                NT const& ratio,
                                NT const& C,
                                NT const& frac,
                                unsigned int const& N,
                                unsigned int const& walk_length,
                                NT const& error,
                                std::vector<NT>& a_vals,
                                MT const& inv_covariance_matrix,
                                RandomNumberGenerator& rng) 
{
    
    typedef typename NonSphericalGaussianFunctor::FunctionFunctor<Point>    Func;
    typedef typename NonSphericalGaussianFunctor::GradientFunctor<Point>    Grad;
    typedef typename NonSphericalGaussianFunctor::HessianFunctor<Point>     Hess;
    typedef typename NonSphericalGaussianFunctor::parameters<NT, Point>     func_params;

    typedef crhmc_input<MT, Point, Func, Grad, ZeroFunctor<Point>> Input;

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

    // compute the first variance
    // for this we need P
    auto ball = P_copy.ComputeInnerBall();
    P_copy.shift(ball.first.getCoefficients()); // when using max_ellipsoid for rounding this center is the origin, but if we use other covariances this is different than the origin
    get_first_gaussian(P_copy, frac, ball.second, error, a_vals); // the function included from volume_cooling_gaussians.hpp (spherical gaussians)
#ifdef VOLESTI_DEBUG
    std::cout << "first gaussian computed " << a_vals[0] << std::endl;
#endif

    //for the rest we need Pin
    NT a_stop = 0.0;
    const NT tol = 0.001;
    unsigned int it = 0;
    unsigned int n = Pin_copy.dimension();
    const unsigned int totalSteps = ((int)150/((1.0 - frac) * error)) + 1;

    if (a_vals[0]<a_stop) a_vals[0] = a_stop;

#ifdef VOLESTI_DEBUG
    std::cout << "Computing the sequence of gaussians..\n" << std::endl;
#endif
    NT initial_eta;
    Point start_point;

    int dimension = Pin_copy.dimension();
    func_params initial_f_params = func_params(Point(dimension), a_vals[0], -1, inv_covariance_matrix);
    Func initial_f(initial_f_params);
    Grad initial_g(initial_f_params);
    ZeroFunctor<Point> initial_zerof;

    Input initial_input = convert2crhmc_input<Input, Polytope, Func, Grad, ZeroFunctor<Point>>(Pin_copy, initial_f, initial_g, initial_zerof);
    CrhmcProblem initial_problem = CrhmcProblem(initial_input);

    Point initial_p = Point(initial_problem.center);
    initial_problem.options.simdLen = simdLen;
    crhmc_walk_params initial_params(initial_input.df, initial_p.dimension(), initial_problem.options);
    CRHMCWalkType initial_walk = CRHMCWalkType(initial_problem, initial_p, initial_input.df, initial_input.f, initial_params);
    
    burnIn<CRHMCWalkType, crhmc_walk_params, simdLen, Grad, Func, CrhmcProblem>(
        initial_p, walk_length, rng, initial_g, initial_f, initial_params, initial_problem, initial_walk, 2.5 * N);

    //fix eta
    initial_eta = initial_walk.get_current_eta();

    //fix point
    start_point = initial_problem.center;

    while (true) {

        // Compute the next gaussian
        NT next_a = get_next_gaussian<simdLen>(
            Pin_copy, start_point, a_vals[it], N, ratio, C, walk_length, rng, inv_covariance_matrix);

        // Decide if this is the last one
        NT curr_fn = 0;
        NT curr_its = 0;
        auto steps = totalSteps;

        //TODO: potential problem creation and preprocessing optimization

        // Create the CRHMC problem for this variance
        int dimension = Pin_copy.dimension();
        func_params f_params = func_params(Point(dimension), a_vals[it], -1, inv_covariance_matrix);

        Func f(f_params);
        Grad g(f_params);
        ZeroFunctor<Point> zerof;

        Input input = convert2crhmc_input<Input, Polytope, Func, Grad, ZeroFunctor<Point>>(Pin_copy, f, g, zerof);

        typedef crhmc_problem<Point, Input> CrhmcProblem;
        
        CrhmcProblem problem = CrhmcProblem(input);

       if(problem.terminate) { return; }

        problem.options.simdLen = simdLen;

        crhmc_walk_params params(input.df, start_point.dimension(), problem.options);

        //eta initialization with the previous walk eta
        problem.options.etaInitialize = false;
        params.eta = initial_eta;

        // Create the walk object for this problem
        CRHMCWalkType walk = CRHMCWalkType(problem, start_point, input.df, input.f, params);

        // Compute some ratios to decide if this is the last gaussian
        for (unsigned int j = 0; j < steps; j++) 
        {
            walk.apply(rng, walk_length);
            start_point = walk.getPoint();
            curr_its += 1.0;    
            curr_fn += eval_exp(start_point, inv_covariance_matrix, next_a, a_vals[it]);
        }

        //restore the new eta and start point, by looking at the walk after its operations
        initial_eta = walk.get_current_eta();  

        // Remove the last gaussian.
        // Set the last a_i equal to zero
        if (next_a > 0 && curr_fn / curr_its > (1.0 + tol)) 
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
}


template <typename NT>
struct non_gaussian_annealing_parameters
{
    non_gaussian_annealing_parameters(unsigned int d)
        : frac(0.1)
        , ratio(NT(1) - NT(1) / NT(d))
        , C(NT(2))
        , N(500 * ((int) C) + ((int) (d * d)))
        , W(6 * d * d + 800)
    {}

    NT frac;
    NT ratio;
    NT C;
    unsigned int N;
    unsigned int W;
};


template
<
    typename Polytope,
    typename RandomNumberGenerator,
    typename WalkTypePolicy = CRHMCWalk,
    int simdLen = 8
>
double non_spherical_crhmc_volume_cooling_gaussians(Polytope& Pin,
                                RandomNumberGenerator& rng,
                                double const& error = 0.1,
                                unsigned int const& walk_length = 1)
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT     NT;
    typedef typename Polytope::VT  VT;
    typedef typename Polytope::MT  MT;
    typedef typename NonSphericalGaussianFunctor::FunctionFunctor<Point>    Func;
    typedef typename NonSphericalGaussianFunctor::GradientFunctor<Point>    Grad;
    typedef typename NonSphericalGaussianFunctor::parameters<NT, Point>     func_params;

    typedef crhmc_input<MT, Point, Func, Grad, ZeroFunctor<Point>> Input;
    typedef crhmc_problem<Point, Input> CrhmcProblem;   

    typedef ImplicitMidpointODESolver<Point, NT, CrhmcProblem, Grad, simdLen> Solver;

    typedef typename WalkTypePolicy::template Walk
            <
                    Point,
                    CrhmcProblem,
                    RandomNumberGenerator,
                    Grad,
                    Func,
                    Solver
            > CRHMCWalkType;

    typedef typename WalkTypePolicy::template parameters
          <
                  NT,
                  Grad
          > crhmc_walk_params;

    typedef HPolytope<Point> HPOLYTOPE;

    HPOLYTOPE P(Pin.dimension(), Pin.get_mat(), Pin.get_vec());
    HPOLYTOPE newPin(Pin.dimension(), Pin.get_mat(), Pin.get_vec());
    unsigned int n = P.dimension();
    unsigned int m = P.num_of_hyperplanes();
    
    //compute inscribed ellipsoid
    NT tol = std::pow(10, -6.0), reg = std::pow(10, -4.0);
    unsigned int maxiter = 100;
    P.normalize();
    VT x0 = compute_feasible_point(P.get_mat(), P.get_vec());
    auto ellipsoid_result = compute_inscribed_ellipsoid<MT, EllipsoidType::MAX_ELLIPSOID>(P.get_mat(), P.get_vec(), x0, maxiter, tol, reg);

    // extract the covariance matrix and the center of the ellipsoid
    MT inv_covariance_matrix = std::get<0>(ellipsoid_result); //this is the covariance to use in the telescopic product
    VT center = std::get<1>(ellipsoid_result);
    MT covariance_matrix = inv_covariance_matrix.inverse(); 

    newPin.shift(center); //we shift the initial polytope so that the origin is the center of the gaussians
    P.shift(center);

    // we apply the rounding transformation
    Eigen::LLT<MT> lltOfA(covariance_matrix);
    auto L = lltOfA.matrixL();
    P.linear_transformIt(L);


    // Initialize the gaussian_annealing_parameters struct
    non_gaussian_annealing_parameters<NT> parameters(P.dimension());

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

    compute_annealing_schedule<simdLen, Point>(newPin, P, ratio, C, parameters.frac, N, walk_length, error, a_vals, inv_covariance_matrix, rng);

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
    
    MT scaled_matrix = a_vals[0] * inv_covariance_matrix;
    NT det_scaled = scaled_matrix.determinant();
    NT vol = std::pow(M_PI, n / 2.0) / std::sqrt(det_scaled);

    unsigned int i=0;

    typedef typename std::vector<NT>::iterator viterator;
    viterator itsIt = its.begin();
    auto avalsIt = a_vals.begin();
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
        int dimension = newPin.dimension();
        func_params f_params = func_params(Point(dimension), *avalsIt, -1, inv_covariance_matrix);

        Func f(f_params);
        Grad g(f_params);

        ZeroFunctor<Point> zerof;

        //create the crhmc problem
        Input input = convert2crhmc_input<Input, Polytope, Func, Grad, ZeroFunctor<Point>>(newPin, f, g, zerof);

        typedef crhmc_problem<Point, Input> CrhmcProblem;
        CrhmcProblem problem = CrhmcProblem(input);
        Point p = problem.center;

        if(problem.terminate){return 0;}
        problem.options.simdLen=simdLen;

        //create the walk and do the burnIn
        crhmc_walk_params params(input.df, p.dimension(), problem.options);
        CRHMCWalkType walk = CRHMCWalkType(problem, p, input.df, input.f, params);

        burnIn<CRHMCWalkType, crhmc_walk_params, simdLen, Grad, Func, CrhmcProblem>(
            p, walk_length, rng, g, f, params, problem, walk, 2.5 * N);

        while (!done || (*itsIt) < min_steps)
        {
            walk.apply(rng, walk_length);
            p = walk.getPoint();
            *itsIt = *itsIt + 1.0;
            
            *fnIt += eval_exp(p, inv_covariance_matrix, *(avalsIt+1), *avalsIt);

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

            if ((max_val - min_val) / max_val <= curr_eps / 2.0)
            {
                done = true;
            }

            index = index % W + 1;
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

#endif // VOLUME_COOLING_NON_SPHERICAL_GAUSSIANS_CRHMC_HPP