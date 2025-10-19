// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef VOLESTI_SIMULATED_ANNEALING_HPP
#define VOLESTI_SIMULATED_ANNEALING_HPP

#include <cmath>
#include <list>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <limits>

#include "generators/boost_random_number_generator.hpp"
#include "random_walks/boltzmann_hmc_walk.hpp"
#include "optimization/sliding_window.hpp"

/// Number of sample points for diameter estimation
/// When estimating the diameter of the spectrahedron,
/// sample 1000 + sqrt(dimension) points to estimate it
#define CONSTANT_1 1000

/// Configuration parameters for the simulated annealing algorithm
/// Contains all tunable parameters for controlling the optimization process,
/// including convergence criteria, temperature schedule, and random walk settings.
/// \tparam Point Point type representing vectors in the optimization space
template<class Point>
struct SimulatedAnnealingSettings {
    /// Numeric type
    typedef typename Point::FT NT;

    /// Desired accuracy threshold for convergence (relative error tolerance)
    NT error;
    
    /// Number of steps taken in each HMC random walk iteration
    int walkLength;
    
    /// Maximum number of optimization steps before termination
    int maxSteps;
    
    /// Temperature decay factor for exponential cooling schedule (must be in (0,1))
    NT decFactor;
    
    /// Ratio of minimum temperature to initial temperature (T_min = T_0 * tempMinRatio)
    NT tempMinRatio;
    
    /// If true, use polynomial cooling schedule; otherwise use exponential schedule
    bool usePolynomialSchedule;
    
    /// Parameter k for polynomial schedule: alpha = 1 - 1/(d*k) where d is dimension
    NT polynomialK;
    
    /// If true, check for early convergence using sliding window analysis
    bool enableEarlyConvergence;
    
    /// Size of sliding window for tracking convergence (number of recent values to track)
    int convergenceWindow;

    /// Construct settings with default or custom parameters
    /// Validates all parameters and throws exceptions if constraints are violated.
    /// \param[in] error_ Convergence tolerance (must be > 0)
    /// \param[in] walkLength_ HMC walk length per iteration (must be > 0)
    /// \param[in] maxSteps_ Maximum optimization steps (must be > 0)
    /// \param[in] decFactor_ Exponential decay factor (must be in (0,1))
    /// \param[in] tempMinRatio_ Minimum temperature ratio (must be > 0)
    /// \param[in] usePolynomialSchedule_ Enable polynomial cooling schedule
    /// \param[in] polynomialK_ Polynomial schedule parameter (must be > 0)
    /// \param[in] enableEarlyConvergence_ Enable early stopping based on convergence
    /// \param[in] convergenceWindow_ Sliding window size for convergence detection
    /// \throws std::invalid_argument if any parameter violates its constraints
    SimulatedAnnealingSettings(NT error_ = NT(1e-6),
                               int walkLength_ = 15,
                               int maxSteps_ = 1000,
                               NT decFactor_ = NT(0.5),
                               NT tempMinRatio_ = NT(1e-8),
                               bool usePolynomialSchedule_ = true,
                               NT polynomialK_ = NT(0.5),
                               bool enableEarlyConvergence_ = true,
                               int convergenceWindow_ = 45)
        : error(error_), walkLength(walkLength_), maxSteps(maxSteps_),
          decFactor(decFactor_), 
          tempMinRatio(tempMinRatio_),
          usePolynomialSchedule(usePolynomialSchedule_),
          polynomialK(polynomialK_),
          enableEarlyConvergence(enableEarlyConvergence_),
          convergenceWindow(convergenceWindow_) {
        if (error <= NT(0))          throw std::invalid_argument("[SimAnn] error must be > 0");
        if (walkLength <= 0)         throw std::invalid_argument("[SimAnn] walkLength must be > 0");
        if (maxSteps <= 0)           throw std::invalid_argument("[SimAnn] maxSteps must be > 0");
        if (decFactor <= NT(0) || decFactor >= NT(1))
            throw std::invalid_argument("[SimAnn] decFactor must be in (0,1)");
        if (tempMinRatio <= NT(0))   throw std::invalid_argument("[SimAnn] tempMinRatio must be > 0");
        if (polynomialK <= NT(0))    throw std::invalid_argument("[SimAnn] polynomialK must be > 0");
    }
};

/// Compute temperature using polynomial cooling schedule
/// Calculates temperature at a given step using the formula:
/// T_i = T_0 * alpha^i where alpha = 1 - 1/(d*k)
/// The decay factor alpha is clamped to [0.5, 0.999] to ensure stable cooling.
/// The computed temperature is always at least T_min.
/// \tparam NT Numeric type for calculations
/// \param[in] step Current optimization step (iteration number)
/// \param[in] T0 Initial temperature
/// \param[in] Tmin Minimum allowed temperature
/// \param[in] dimension Dimension of the optimization space
/// \param[in] k Polynomial schedule parameter controlling cooling rate
/// \return Temperature at the given step, constrained to be >= Tmin
template <typename NT>
NT computePolynomialTemperature(int step, NT T0, NT Tmin, int dimension, NT k) {
    // Compute decay factor: alpha = 1 - 1/(d*k)
    NT alpha = NT(1) - NT(1) / (NT(dimension) * k);
    
    // Clamp alpha to safe range to prevent too fast or too slow cooling
    alpha = std::max(NT(0.5), std::min(alpha, NT(0.999)));
    
    // Apply polynomial schedule: T = T_0 * alpha^step
    NT T = T0 * std::pow(alpha, NT(step));
    
    // Ensure temperature doesn't fall below minimum
    return std::max(T, Tmin);
}

/// Solve semidefinite programming problem using simulated annealing
/// Minimizes a linear objective function c^T * x subject to a linear matrix
/// inequality constraint (spectrahedron). Uses Hamiltonian Monte Carlo (HMC)
/// random walk with Boltzmann distribution for sampling, combined with simulated
/// annealing temperature schedule for optimization.
/// Algorithm outline:
/// 1. Normalize objective function and estimate feasible region diameter
/// 2. Initialize temperature schedule (T_0 = diameter)
/// 3. At each step:
///    - Update temperature according to schedule (polynomial or exponential)
///    - Sample new point using HMC with current temperature
///    - Update best solution if improvement found
///    - Check for convergence if early stopping enabled
/// 4. Return best solution found
/// \tparam Spectrahedron Type representing the feasible region (linear matrix inequality)
/// \tparam Point Point type representing vectors in optimization space
/// \tparam Settings Settings type (must be SimulatedAnnealingSettings or compatible)
/// \param[in] spectrahedron The feasible region defined by LMI constraints
/// \param[in] objectiveFunction Linear objective function to minimize (c^T * x)
/// \param[in] settings Algorithm parameters controlling optimization behavior
/// \param[in] interiorPoint Initial feasible solution (must be strictly interior)
/// \param[out] solution Output parameter storing the best solution found
/// \param[in] verbose If true, print progress information during optimization
/// \return Objective function value at the best solution found
/// \throws std::runtime_error if objective function has zero norm or diameter estimation fails
template <typename Spectrahedron, typename Point, typename Settings>
typename Point::FT solve_sdp(Spectrahedron& spectrahedron,
                             Point const& objectiveFunction,
                             Settings const& settings,
                             Point const& interiorPoint,
                             Point& solution,
                             bool verbose = false) {

    typedef typename Spectrahedron::NT NT;
    typedef typename Spectrahedron::VT VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef BoltzmannHMCWalk::Walk<Spectrahedron, RNGType> HMC;

    // Extract and normalize objective function coefficients
    VT c_raw = objectiveFunction.getCoefficients();
    NT cn_raw = c_raw.norm();
    if (cn_raw <= NT(0)) 
        throw std::runtime_error("[SimAnn] Objective has zero norm");

    // Create normalized version for HMC (prevents numerical issues)
    VT c_unit = c_raw / cn_raw;
    Point objective_for_hmc(c_unit);

    // Initialize random number generator
    RNGType rng(spectrahedron.dimension());

    // Estimate diameter of feasible region by sampling
    // Number of samples: CONSTANT_1 + sqrt(dimension)
    NT diameter = spectrahedron.estimateDiameter(
        CONSTANT_1 + std::sqrt(spectrahedron.dimension()),
        interiorPoint, 
        rng
    );

    if (diameter <= NT(0)) 
        throw std::runtime_error("[SimAnn] Invalid diameter estimation");

    // Caution: Empirical modifications about diameter and initial temperature
    const unsigned int dim = spectrahedron.dimension();

    if(dim >= 100) 
        diameter = std::max(diameter, NT(1e10));
    else 
        diameter = std::max(diameter, NT(1));

    // Initial temperature set to estimated diameter
    NT T0 = std::min(diameter, NT(1.0));


    // Minimum temperature as fraction of initial temperature
    NT Tmin = T0 * settings.tempMinRatio;

    if (verbose) {
        std::cout << "[SimAnn] T0=" << T0 << ", Tmin=" << Tmin 
                  << ", diameter=" << diameter << std::endl;
        std::cout << "[SimAnn] Cooling: " 
                  << (settings.usePolynomialSchedule ? "polynomial" : "exponential");
        if (settings.usePolynomialSchedule) {
            NT alpha = NT(1) - NT(1) / (NT(spectrahedron.dimension()) * settings.polynomialK);
            alpha = std::min(alpha, NT(0.999));
            std::cout << " (k=" << settings.polynomialK << ", alpha=" << alpha << ")";
        } else {
            std::cout << " (alpha=" << settings.decFactor << ")";
        }
        std::cout << std::endl;
    }

    // Configure Hamiltonian Monte Carlo sampler
    typename HMC::Settings hmc_settings(
        settings.walkLength,
        rng,
        objective_for_hmc,
        T0,
        diameter,
        1000,
        NT(1.0)
    );

    HMC hmcRandomWalk(hmc_settings);

    // State initialization
    // Current point in optimization
    Point current = interiorPoint;
    // Best point found so far
    Point best = current;

    // Evaluate objective at initial point (using unnormalized coefficients)
    NT fCurrent = c_raw.dot(current.getCoefficients());
    NT fBest = fCurrent;

    if (verbose) {
        std::cout << "[SimAnn] Initial objective: " << fBest << std::endl;
    }

    // Initialize sliding window for convergence tracking
    SlidingWindow<NT> convergenceWindow(settings.convergenceWindow);

    // Buffer for storing sampled points
    std::list<Point> buf;
    NT T = T0;

    //Main optimization loop
    for (int step = 0; step < settings.maxSteps; ++step) {
        
        // Update temperature according to selected cooling schedule
        if (settings.usePolynomialSchedule) {
            T = computePolynomialTemperature(step, T0, Tmin,
                                            spectrahedron.dimension(), 
                                            settings.polynomialK);
        } else {
            // Exponential cooling: T_new = T_old * decFactor
            if (step > 0) T *= settings.decFactor;
        }
        
        // Check if minimum temperature reached
        if (T < Tmin) {
            if (verbose) {
                std::cout << "[SimAnn] Minimum temperature " << Tmin 
                          << " reached at step " << step << std::endl;
            }
            break;
        }
        
        // Sample new point using HMC with current temperature
        hmcRandomWalk.setTemperature(T);
        buf.clear();
        hmcRandomWalk.apply(spectrahedron, current, 1, buf);

        if (!buf.empty()) {
            // Update current point with sampled point
            current = buf.back();
            fCurrent = c_raw.dot(current.getCoefficients());
            
            // Update best solution if significant improvement found
            if (fCurrent < fBest - settings.error) {
                best = current;
                fBest = fCurrent;
                if (verbose) {
                    std::cout << "[SimAnn] step " << step
                              << " | best=" << fBest 
                              << " | T=" << T << std::endl;
                }
            }
            
            // Track best value in sliding window for convergence detection
            if (settings.enableEarlyConvergence) {
                convergenceWindow.push(fBest);
            }
        }
        
        // Check for early convergence using sliding window statistics
        if (settings.enableEarlyConvergence && convergenceWindow.isFull()) {
            NT relError = convergenceWindow.getRelativeError();
            
            // Stop if relative change in objective is below threshold
            if (relError < settings.error) {
                if (verbose) {
                    std::cout << "[SimAnn] Early convergence at step " << step 
                              << " (relative error: " << relError << ")" << std::endl;
                }
                break;
            }
        }
    }

    // Final output
    if (verbose) {
        NT final_acceptance = hmcRandomWalk.getAcceptanceRate();
        std::cout << "[SimAnn] Optimization completed" << std::endl;
        std::cout << "[SimAnn] Final objective: " << fBest << std::endl;
        std::cout << "[SimAnn] Final temperature: " << T << std::endl;
        std::cout << "[SimAnn] HMC acceptance: " 
                  << static_cast<double>(final_acceptance * 100.0) << "%" << std::endl;
    }

    // Store best solution in output parameter
    solution = best;
    return fBest;
}

#endif // VOLESTI_SIMULATED_ANNEALING_HPP