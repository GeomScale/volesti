// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Angelos Korakitis, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SIMULATED_ANNEALING_HPP
#define VOLESTI_SIMULATED_ANNEALING_HPP

#include <cmath>
#include <list>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "generators/boost_random_number_generator.hpp"
#include "random_walks/boltzmann_hmc_walk.hpp"

/// A magic number!
/// When estimating the diameter of the spectrahedron,
/// sample 2000 + sqrt(dimension) points to estimate it
#define CONSTANT_1 2000

/// Holds parameters of the simulated annealing algorithm
template<class Point>
struct SimulatedAnnealingSettings {
    /// The numeric type
    typedef typename Point::FT NT;

    /// Desired accuracy (tolerance on improvement - used as soft stop)
    NT error;
    /// The walk length of the HMC random walk (burn steps per sample)
    int walkLength;
    /// Maximum number of SA iterations
    int maxSteps;
    /// Exponential cooling factor (must be in (0,1))
    NT decFactor;
    /// Minimum temperature ratio: Tmin = T0 * tempMinRatio
    NT tempMinRatio;

    SimulatedAnnealingSettings(NT error_ = NT(1e-6),
                               int walkLength_ = 10,
                               int maxSteps_ = 1000,
                               NT decFactor_ = NT(0.99),
                               NT tempMinRatio_ = NT(1e-12))
        : error(error_), walkLength(walkLength_), maxSteps(maxSteps_),
          decFactor(decFactor_), tempMinRatio(tempMinRatio_) {
        if (error <= NT(0))          throw std::invalid_argument("error must be > 0");
        if (walkLength <= 0)         throw std::invalid_argument("walkLength must be > 0");
        if (maxSteps <= 0)           throw std::invalid_argument("maxSteps must be > 0");
        if (decFactor <= NT(0) || decFactor >= NT(1))
            throw std::invalid_argument("decFactor must be in (0,1)");
        if (tempMinRatio <= NT(0))   throw std::invalid_argument("tempMinRatio must be > 0");
    }
};

/// Helper function for feasibility checking
/// \param[in] S The spectrahedron
/// \param[in] p The point to check
/// \return True if the point is feasible (interior or on boundary), false otherwise
template <typename Spectrahedron, typename Point>
static inline bool is_feasible(const Spectrahedron& S, const Point& p) {
    try { 
        return !S.isExterior(p.getCoefficients()); 
    }
    catch (...) { 
        return false; 
    }
}

/// Simulated Annealing algorithm for semidefinite programming
/// Minimizes c^T x subject to LMI(x) >= 0
/// \param[in] spectrahedron A spectrahedron described by a linear matrix inequality
/// \param[in] objectiveFunction The objective function vector c to minimize
/// \param[in] settings Parameters of the algorithm
/// \param[in] interiorPoint An initial feasible solution to start the algorithm
/// \param[out] solution The vector minimizing the objective function
/// \param[in] verbose True to print diagnostic messages. Default is false
/// \return The optimal objective value found
template <typename Spectrahedron, typename Point, typename Settings>
typename Point::FT solve_sdp(Spectrahedron& spectrahedron,
                             Point const& objectiveFunction,
                             Settings const& settings,
                             Point const& interiorPoint,
                             Point& solution,
                             bool verbose = false) {

    // Fetch the data types we will use
    typedef typename Spectrahedron::NT NT;
    typedef typename Spectrahedron::VT VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef BoltzmannHMCWalk::Walk<Spectrahedron, RNGType> HMC;

    // Validate initial feasibility
    if (!is_feasible(spectrahedron, interiorPoint))
        throw std::runtime_error("Initial point infeasible");

    /******** Objective function handling *********/
    // Keep RAW objective for reporting (fixes scale mismatch vs. test harness)
    // Use UNIT objective for HMC dynamics to ensure proper temperature scaling
    VT c_raw = objectiveFunction.getCoefficients();
    NT cn_raw = c_raw.norm();
    if (cn_raw <= NT(0)) 
        throw std::runtime_error("Objective has zero norm");

    VT c_unit = c_raw / cn_raw;                 // Normalized: ||c_unit|| = 1
    Point objective_for_hmc(c_unit);            // Unit objective for HMC dynamics

    // Initialize random number generator
    RNGType rng(spectrahedron.dimension());

    /******** Diameter estimation *********/
    // Estimate the diameter of the spectrahedron
    // Needed for the random walk and temperature initialization
    NT diameter = spectrahedron.estimateDiameter(
        CONSTANT_1 + std::sqrt(spectrahedron.dimension()),
        interiorPoint, 
        rng
    );
    if (diameter <= NT(0)) 
        throw std::runtime_error("Invalid diameter estimation");

    /******** Temperature initialization *********/
    // With unit objective inside HMC, T0 ~ diameter is a natural energy scale
    NT T0   = std::max(NT(1), diameter);
    NT Tmin = T0 * settings.tempMinRatio;

    if (verbose) {
        std::cout << "[SimAnn] T0=" << T0 
                  << "  Tmin=" << Tmin 
                  << "  diameter=" << diameter << std::endl;
    }

    /******** Initialize HMC random walk *********/
    // HMC settings with unit objective for proper energy scaling
    typename HMC::Settings hmc_settings(
        settings.walkLength,           // walk length (burn-in steps)
        rng,                          // random number generator
        objective_for_hmc,            // unit-normalized objective
        T0,                           // initial temperature
        diameter,                     // estimated diameter
        10000,                        // reflections cap (increased for stability)
        NT(1.0)                       // trajectory_factor
    );

    HMC hmcRandomWalk(hmc_settings);

    /******** State initialization *********/
    Point current = interiorPoint;
    Point best    = current;

    // Evaluate with RAW objective for correctness vs. external scales
    NT fCurrent  = c_raw.dot(current.getCoefficients());
    NT fBest     = fCurrent;

    if (verbose) {
        std::cout << "[SimAnn] Initial objective value: " << fBest << std::endl;
    }

    std::list<Point> buf;
    NT T = T0;

    /******** Main optimization loop *********/
    for (int step = 0; step < settings.maxSteps; ++step) {
        
        // Sample one point at temperature T via HMC
        hmcRandomWalk.setTemperature(T);
        buf.clear();
        hmcRandomWalk.apply(spectrahedron, current, 1, buf);

        if (!buf.empty()) {
            current = buf.back();

            // Evaluate with RAW objective (kept consistent for reporting/return)
            fCurrent = c_raw.dot(current.getCoefficients());
            
            // Update best solution if improved beyond tolerance
            if (fCurrent < fBest - settings.error) {
                best = current;
                fBest = fCurrent;
                if (verbose) {
                    std::cout << "[SimAnn] step " << step
                              << " | new best=" << fBest 
                              << " | T=" << T << std::endl;
                }
            }
        }

        // Apply exponential cooling schedule
        T *= settings.decFactor;
        
        // Check minimum temperature termination condition
        if (T < Tmin) {
            if (verbose) {
                std::cout << "[SimAnn] Minimum temperature " << Tmin 
                          << " reached at step " << step << std::endl;
            }
            break;
        }
    }

    /******** Final output *********/
    if (verbose) {
        NT final_hmc_acceptance = hmcRandomWalk.getAcceptanceRate();
        std::cout << "[SimAnn] Optimization completed" << std::endl;
        std::cout << "[SimAnn] Final objective value: " << fBest << std::endl;
        std::cout << "[SimAnn] Final temperature: " << T << std::endl;
        std::cout << "[SimAnn] HMC acceptance rate: " 
                  << static_cast<double>(final_hmc_acceptance * 100.0) << "%" << std::endl;
    }

    solution = best;
    // Return RAW objective value (not normalized)
    return fBest;
}

/// Convenience overload with default settings
/// \param[in] spectrahedron A spectrahedron described by a linear matrix inequality
/// \param[in] objectiveFunction The objective function vector c to minimize
/// \param[in] interiorPoint An initial feasible solution to start the algorithm
/// \param[out] solution The vector minimizing the objective function
/// \param[in] verbose True to print diagnostic messages. Default is false
/// \return The optimal objective value found
template <typename Spectrahedron, typename Point>
inline typename Point::FT solve_sdp(Spectrahedron& spectrahedron,
                                   Point const& objectiveFunction,
                                   Point const& interiorPoint,
                                   Point& solution,
                                   bool verbose = false) {
    SimulatedAnnealingSettings<Point> defaultSettings;
    return solve_sdp(spectrahedron, objectiveFunction, defaultSettings, 
                     interiorPoint, solution, verbose);
}

#endif // VOLESTI_SIMULATED_ANNEALING_HPP