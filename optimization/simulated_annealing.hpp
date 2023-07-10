// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SIMULATED_ANNEALING_HPP
#define VOLESTI_SIMULATED_ANNEALING_HPP

#include <format.hpp>

#include "optimization/sliding_window.hpp"
#include "random_walks/boltzmann_hmc_walk.hpp"


/// A magic number!
/// when estimating the diameter of the spectrahedron,
/// sample 20 + sqrt(dimension) points to estimate it
#define CONSTANT_1 20

/// Holds parameters of the algorithm
/// \tparam Point Point Type
template<class Point>
struct SimulatedAnnealingSettings {
    /// The numeric type
    typedef typename Point::FT NT;

    /// Desired accuracy (relative error)
    NT error;
    /// The walk length of the HMC random walk
    int walkLength;
    /// A bound to the number of steps; if negative it is unbounded
    int maxNumSteps;

    /// Starting from an initial temperature, at each step it will decrease by a factor of
    /// \[ 1 - 1 / (dimension^k) \]. Default is 0.5
    NT k;

    SimulatedAnnealingSettings(NT const error, int const walkLength = 1, int const maxNumSteps = -1, NT const k = 0.5) : error(error),
        walkLength(walkLength), maxNumSteps(maxNumSteps), k(k) {}
};


/// Simulated Annealing algorithm for a semidefinite program
/// Minimize \[ c^T x \], s.t. LMI(x) <= 0
/// \param[in] spectrahedron A spectrahedron described by a linear matrix inequality
/// \param[in] objectiveFunction The function we minimize
/// \param[in] settings Parameters of the algorithm
/// \param[in] interiorPoint An initial feasible solution to start the algorithm
/// \param[out] solution The vector minimizing the objective function
/// \param[in] verbose True to print messages. Default is false
/// \return The best approximation to the optimal solution
template <typename _Spectrahedron, typename Point, typename _Settings>
double solve_sdp(_Spectrahedron & spectrahedron, Point const & objectiveFunction, _Settings const & settings,
         Point const & interiorPoint, Point& solution, bool verbose = false) {

    // fetch the data types we will use
    typedef  typename _Spectrahedron::NT NT;
    typedef  typename _Spectrahedron::MT MT;
    typedef  typename _Spectrahedron::VT VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef BoltzmannHMCWalk::Walk<_Spectrahedron, RNGType > HMC;

    // the algorithm requires the objective function to be normalized
    // we will need to remember the norm
    VT _objectiveFunctionNormed = objectiveFunction.getCoefficients();
    NT objectiveFunctionNorm = _objectiveFunctionNormed.norm();
    _objectiveFunctionNormed.normalize();
    Point objectiveFunctionNormed = Point(_objectiveFunctionNormed);

    RNGType rng(spectrahedron.dimension());

    // Estimate the diameter of the spectrahedron
    // needed for the random walk and for the simulated annealing algorithm
    NT diameter = spectrahedron.estimateDiameter(CONSTANT_1 + std::sqrt(spectrahedron.dimension()), interiorPoint, rng);

    /******** initialization *********/
    solution = interiorPoint;
    // the minimum till last iteration
    NT currentMin = objectiveFunction.dot(solution);
    int stepsCount = 0;
    // initial temperature must be the diameter of the body
    NT temperature = diameter;
    // after each iteration, temperature = temperature * tempDecreaseFactor
    NT tempDecreaseFactor = 1.0 - static_cast<NT>(1.0 / std::pow(spectrahedron.dimension(), settings.k));

    // initialize random walk;
    typename HMC::Settings hmc_settings = typename HMC::Settings(settings.walkLength, rng, objectiveFunction, temperature, diameter);
    HMC hmcRandomWalk = HMC(hmc_settings);
    NT previous_min = objectiveFunction.dot(solution);

    /******** solve *********/
    // if settings.maxNumSteps is negative there is no
    // bound to the number of steps - stop
    // when desired relative error is achieved
    while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) {

        // sample one point with current temperature
        std::list<Point> randPoints;

        // get a sample under the Boltzmann distribution
        // using the HMC random walk
        while (1) {
            hmcRandomWalk.apply(spectrahedron, solution, settings.walkLength, randPoints);

            // if the sampled point is not inside the spectrahedron (error in boundary oracle),
            // get a new one
            if (spectrahedron.isExterior(spectrahedron.get_C())) {
                if (verbose) std::cout << "Sampled point outside the spectrahedron.\n";
                randPoints.clear();
                spectrahedron.resetFlags();
            }
            else {
                // update values;
                solution = randPoints.back();
                randPoints.clear();
                break;
            }
        }

        // update current value
        currentMin = objectiveFunction.dot(solution);
        ++stepsCount;

        // compute relative error
        NT relError = relativeError(previous_min, currentMin);
        previous_min = currentMin;

        if (verbose)
            std::cout << "Step: " << stepsCount << ", Temperature: " << temperature << ", Min: " << currentMin
                      << ", Relative error: " << relError << "\n";

        // check if we reached desired accuracy
        if (relError < settings.error)
            break;

        // decrease the temperature
        temperature *= tempDecreaseFactor;
        hmcRandomWalk.setTemperature(temperature);

    } /* while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) { */

    // return the minimum w.r.t. the original objective function
    return currentMin*objectiveFunctionNorm;
}



#endif //VOLESTI_SIMULATED_ANNEALING_HPP
