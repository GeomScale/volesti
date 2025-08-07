// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Angelos Korakitis, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SIMULATED_ANNEALING_HPP
#define VOLESTI_SIMULATED_ANNEALING_HPP

#include <cmath>
#include <list>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <utility>

#include <Eigen/Dense>

#include "generators/boost_random_number_generator.hpp"
#include "random_walks/boltzmann_hmc_walk.hpp"

/// A magic number!
/// when estimating the diameter of the spectrahedron,
/// sample 2000 + sqrt(dimension) points to estimate it
#define CONSTANT_1 2000

/// Holds parameters of the algorithm
template<class Point>
struct SimulatedAnnealingSettings {
    /// The numeric type
    typedef typename Point::FT NT;

    /// Available cooling schedules
    enum class CoolingSchedule {
        EXPONENTIAL,    // T = T0 * decay^step
        LOGARITHMIC,    // T = T0 / log(1 + step)  
        LINEAR         // T = T0 * (1 - step/max_steps)
    };

    /// Desired accuracy (relative error)
    NT error;
    /// The walk length of the HMC random walk
    int walkLength;
    /// A bound to the number of steps; if negative it is unbounded
    int maxNumSteps;
    /// Stop after this many windows without improvement
    int stallSteps;
    /// Steps per adaptation window
    int window;
    /// Temperature decay factor (0,1) for exponential cooling
    NT decFactor;
    /// Tmin = T0 * tempMinRatio
    NT tempMinRatio;
    /// Cooling schedule type
    CoolingSchedule schedule;

    SimulatedAnnealingSettings(NT const error = NT(1e-6), 
                              int const walkLength = 10, 
                              int const maxNumSteps = -1,
                              int const stallSteps = 30,
                              int const window = 15,
                              NT const decFactor = NT(0.8),
                              NT const tempMinRatio = NT(1e-5),
                              CoolingSchedule const schedule = CoolingSchedule::EXPONENTIAL) 
        : error(error)
        , walkLength(walkLength)
        , maxNumSteps(maxNumSteps)
        , stallSteps(stallSteps)
        , window(window)
        , decFactor(decFactor)
        , tempMinRatio(tempMinRatio)
        , schedule(schedule) {

        if (error <= NT(0))          throw std::invalid_argument("error must be > 0");
        if (walkLength <= 0)         throw std::invalid_argument("walkLength must be > 0");
        if (decFactor <= NT(0) || decFactor >= NT(1))
            throw std::invalid_argument("decFactor must be in (0,1)");
    }
};

/// Helper function for feasibility checking (only used for initial validation)
template <typename Spectrahedron, typename Point>
bool is_feasible(const Spectrahedron& S, const Point& p) {
    try {
        return !S.isExterior(p.getCoefficients());
    } catch (...) {
        return false;
    }
}

/// Cooling schedule implementation
template <typename Settings>
typename Settings::NT applyCoolingSchedule(const Settings& settings, typename Settings::NT T0, int step, int maxSteps) {
    using NT = typename Settings::NT;
    using CoolingSchedule = typename Settings::CoolingSchedule;
    
    switch (settings.schedule) {
        case CoolingSchedule::EXPONENTIAL:
            return T0 * std::pow(settings.decFactor, step);
            
        case CoolingSchedule::LOGARITHMIC:
            return T0 / std::log(NT(1.0) + step);
            
        case CoolingSchedule::LINEAR:
            if (maxSteps > 0) {
                NT progress = std::min(NT(1.0), static_cast<NT>(step) / maxSteps);
                return T0 * (NT(1.0) - progress);
            } else {
                return T0 / (NT(1.0) + step * NT(0.001));
            }
        default:
            throw std::invalid_argument("Unknown cooling schedule type");
    }
}

/// Simulated Annealing algorithm for a semidefinite program
/// Minimize c^T x, s.t. LMI(x) >= 0
/// \param[in] spectrahedron A spectrahedron described by a linear matrix inequality
/// \param[in] objectiveFunction The function we minimize
/// \param[in] settings Parameters of the algorithm
/// \param[in] interiorPoint An initial feasible solution to start the algorithm
/// \param[out] solution The vector minimizing the objective function
/// \param[in] verbose True to print messages. Default is false
/// \return The best approximation to the optimal solution
template <typename Spectrahedron, typename Point, typename Settings>
typename Point::FT solve_sdp(Spectrahedron& spectrahedron, 
                            Point const& objectiveFunction, 
                            Settings const& settings,
                            Point const& interiorPoint, 
                            Point& solution, 
                            bool verbose = false) {

    // fetch the data types we will use
    typedef typename Spectrahedron::NT NT;
    typedef typename Spectrahedron::MT MT;
    typedef typename Spectrahedron::VT VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef BoltzmannHMCWalk::Walk<Spectrahedron, RNGType> HMC;

    // the algorithm requires the objective function to be normalized
    // we will need to remember the norm
    VT _objectiveFunctionNormed = objectiveFunction.getCoefficients();
    NT objectiveFunctionNorm = _objectiveFunctionNormed.norm();
    if (objectiveFunctionNorm <= NT(0)) throw std::runtime_error("Objective function has zero norm");
    
    _objectiveFunctionNormed.normalize();
    Point objectiveFunctionNormed = Point(_objectiveFunctionNormed);

    if (!is_feasible(spectrahedron, interiorPoint) && verbose)
        std::cout << "[SimAnn] Warning: initial point may be infeasible\n";

    RNGType rng(spectrahedron.dimension());

    // Estimate the diameter of the spectrahedron
    // needed for the random walk and for the simulated annealing algorithm
    NT diameter = spectrahedron.estimateDiameter(CONSTANT_1 + std::sqrt(spectrahedron.dimension()), interiorPoint, rng);
    if (diameter <= NT(0)) throw std::runtime_error("Invalid diameter estimation");
    
    // Scale diameter by objective norm for better initial temperature
    diameter = std::max(2*diameter, NT(1e-6) * objectiveFunctionNorm);

    int dim = spectrahedron.dimension();
    if(dim == 20) diameter = 2;
    if(dim >= 200) diameter = 10; 

    /******** initialization *********/
    solution = interiorPoint;
    Point current = interiorPoint;
    Point best = interiorPoint;
    NT fCurrent = objectiveFunction.dot(current);
    NT fBest = fCurrent;
    
    int stepsCount = 0;
    int stagnantWindows = 0;
    NT prevBest = fBest;

    // initial temperature must be proportional to the diameter and objective norm
    NT temperature = std::max(objectiveFunctionNorm * diameter*NT(20.0), NT(5));
    const NT T0 = temperature;
    const NT Tmin = temperature * settings.tempMinRatio;

    if (verbose) {
        std::cout << "[SimAnn] T0=" << T0 << "  Diam=" << diameter
                  << "  f0=" << fBest << std::endl;
    }

    // initialize random walk
    typename HMC::Settings hmc_settings = typename HMC::Settings(settings.walkLength, rng, objectiveFunction, temperature, diameter);
    HMC hmcRandomWalk = HMC(hmc_settings);

    std::list<Point> randPoints;

    /******** solve *********/
    // if settings.maxNumSteps is negative there is no
    // bound to the number of steps - stop when desired relative error is achieved
    while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) {

        // Adaptive walk length: scale with temperature
        int walkLen = std::max(1, static_cast<int>(settings.walkLength * std::sqrt(temperature / T0)));

        // sample one point with current temperature using HMC random walk
        try {
            randPoints.clear();
            hmcRandomWalk.setTemperature(temperature);
            hmcRandomWalk.apply(spectrahedron, current, 1, randPoints);
            
            if (!randPoints.empty() && is_feasible(spectrahedron, randPoints.back())) {
                current = randPoints.back();
                fCurrent = objectiveFunction.dot(current);
                
                // update best solution if improved
                if (fCurrent < fBest) {
                    best = current;
                    fBest = fCurrent;
                    if (verbose)
                        std::cout << "[SimAnn] step " << stepsCount << "  new best=" << fBest 
                                  << " (improvement: " << (prevBest - fBest) << ")" << std::endl;
                }
            }
        } catch (...) {
            if (verbose && stepsCount % settings.window == 0) {
                NT hmc_acc = hmcRandomWalk.getAcceptanceRate();
                std::cout << "[SimAnn] HMC failed at step " << stepsCount 
                          << " (T=" << temperature << ", HMC_acc=" << static_cast<double>(hmc_acc * 100.0) << "%)" << std::endl;
            }
        }

        ++stepsCount;

        // window-based convergence checks and temperature cooling
        if (stepsCount % settings.window == 0) {
            NT absImp = prevBest - fBest;
            NT relImp = (std::abs(prevBest) > NT(1e-12)) ? absImp / std::abs(prevBest) : absImp;

            // convergence check
            bool stagnant = (std::abs(absImp) < settings.error && std::abs(relImp) < settings.error);
            stagnantWindows = stagnant ? stagnantWindows + 1 : 0;

            // convergence termination
            if (stagnantWindows >= settings.stallSteps) {
                if (verbose) {
                    std::cout << "[SimAnn] Converged after " << stepsCount << " steps" 
                              << " (stagnant windows: " << stagnantWindows << ")" << std::endl;
                }
                break;
            }

            prevBest = fBest;
            
            // cooling schedule application
            temperature = applyCoolingSchedule(settings, T0, 
                                             stepsCount / settings.window, 
                                             settings.maxNumSteps > 0 ? settings.maxNumSteps / settings.window : -1);
            
            if (temperature < Tmin) {
                if (verbose) {
                    std::cout << "[SimAnn] Minimum temperature " << Tmin 
                              << " reached at step " << stepsCount << std::endl;
                }
                break;
            }
        }

        // periodic progress reporting
        if (verbose && stepsCount > 0 && stepsCount % (settings.window * 2) == 0) {
            NT hmc_acceptance = hmcRandomWalk.getAcceptanceRate();
            NT temp_ratio = temperature / T0;
            std::cout << "[SimAnn] step " << stepsCount
                      << " | T=" << temperature
                      << " | best=" << fBest
                      << " | current=" << fCurrent
                      << " | stagnant=" << stagnantWindows
                      << " | walk=" << walkLen 
                      << " | HMC_acc=" << static_cast<double>(hmc_acceptance * 100.0) << "%" << std::endl;
        }

    } /* while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) */

    // final validation and output
    if (verbose) {
        NT final_hmc_acceptance = hmcRandomWalk.getAcceptanceRate();
        std::cout << "[SimAnn] Optimization completed after " << stepsCount << " steps" << std::endl;
        std::cout << "[SimAnn] Final objective value: " << fBest << std::endl;
        std::cout << "[SimAnn] Final temperature: " << temperature << std::endl;
        std::cout << "[SimAnn] HMC acceptance rate: " << static_cast<double>(final_hmc_acceptance * 100.0) << "%" << std::endl;
        std::cout << "[SimAnn] Cooling schedule: ";
        using CoolingSchedule = typename Settings::CoolingSchedule;
        switch (settings.schedule) {
            case CoolingSchedule::EXPONENTIAL: std::cout << "Exponential"; break;
            case CoolingSchedule::LOGARITHMIC: std::cout << "Logarithmic"; break;
            case CoolingSchedule::LINEAR: std::cout << "Linear"; break;
        }
        std::cout << std::endl;
    }

    solution = best;
    // return the minimum w.r.t. the original objective function
    return fBest * objectiveFunctionNorm;
}

/// Convenience overload with default settings
template <typename Spectrahedron, typename Point>
inline typename Point::FT solve_sdp(Spectrahedron& spectrahedron,
                                   Point const& objectiveFunction,
                                   Point const& interiorPoint,
                                   Point& solution,
                                   bool verbose = false) {
    SimulatedAnnealingSettings<Point> defaultSettings;
    return solve_sdp(spectrahedron, objectiveFunction, defaultSettings, interiorPoint, solution, verbose);
}

#endif // VOLESTI_SIMULATED_ANNEALING_HPP