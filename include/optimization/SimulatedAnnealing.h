//
// Created by panagiotis on 2/28/20.
//

#ifndef VOLESTI_SIMULATEDANNEALING_H
#define VOLESTI_SIMULATEDANNEALING_H

#include "SlidingWindow.h"

/// Simulated Annealing algorithm for the semidefinite program
/// Minimize \[ c^T x \], s.t. LMI(x) <= 0
/// \tparam Point Class point
/// \tparam MT Matrix type
/// \tparam VT Vector type
template<class Point, typename MT, typename VT>
class SimulatedAnnealing {
public:

    /// The numeric type
    typedef typename Point::FT NT;
    /// The type of the spectrahedron
    typedef Spectrahedron<NT, MT, VT> SPECTRAHEDRON;
    /// To use to generate random numbers
    typedef boost::mt19937 RNGType;
    /// The Hamiltonian Monte Carlo random walk with reflections
    typedef HMC_RandomWalk<Point, MT, VT, RNGType > HMC;

    /// Holds parameters of the algorithm
    struct Settings {
        /// Desired accuracy (relative error)
        NT error;
        /// The walk length of the HMC random walk
        int walkLength;
        /// A bound to the number of steps; if negative it is unbounded
        int maxNumSteps;

        /// Starting from an initial temperature, at each step it will decrease by a factor of
        /// \[ 1 - 1 / (dimension^k) \]. Default is 0.5
        NT k;

        Settings(NT const error, int const walkLength, int const maxNumSteps = -1, NT const k = 0.5) : error(error),
            walkLength(walkLength), maxNumSteps(maxNumSteps), k(k) {}
    };


    /// The feasible region of the semidefinite program
    SPECTRAHEDRON *spectrahedron;
    /// The normalized objective function of the semidefinite program
    Point objectiveFunction;
    /// The norm of the initial objective function
    NT objectiveFunctionNorm;
    /// The parameters of the algorithm
    Settings settings;
    /// A feasible point to start the algorithm
    Point* interiorPoint;
    /// True if we computed the interior point
    bool computedInteriorPoint;
    /// The diameter of the spectrahedron
    NT diameter;

    SimulatedAnnealing() {}

    /// Construct and initialize an instance of this class
    /// \param[in] spectrahedron A spectrahedron described by a linear matrix inequality
    /// \param[in] objectiveFunction The function we minimize
    /// \param[in] settings Parameters of the algorithm
    /// \param[in] interiorPoint An initial feasible solution to start the algorithm
    SimulatedAnnealing(SPECTRAHEDRON* spectrahedron, Point const & objectiveFunction, Settings const & settings,
            Point* interiorPoint = nullptr) : spectrahedron(spectrahedron), settings(settings),
            interiorPoint(interiorPoint) {

        // check if we must compute an interior point
        if (interiorPoint != nullptr)
            computedInteriorPoint = false;
        else {
            this->interiorPoint = new Point(spectrahedron->dimension());
            computedInteriorPoint = true;
        }

        // the algorithm requires the objective function to be normalized
        VT _objectiveFunction = objectiveFunction.getCoefficients();
        objectiveFunctionNorm = _objectiveFunction.norm();
        _objectiveFunction.normalize();
        this->objectiveFunction = Point(_objectiveFunction);

        // Estimate the diameter of the spectrahedron
        // needed for the random walk
        diameter = spectrahedron->estimateDiameter(20 + std::sqrt(spectrahedron->dimension()), *(this->interiorPoint));
    }


    /// Initialize the hamiltonian monte carlo random walk
    /// \param[out] hmcRandomWalk
    /// \param[in] temperature
    void initializeHMC(HMC& hmcRandomWalk, NT const temperature) {
        typedef boost::mt19937 RNGType;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);

        typename HMC::Settings settings = typename HMC::Settings(this->settings.walkLength, rng, objectiveFunction, temperature, diameter);
        hmcRandomWalk = HMC(settings);
    }

    /// Solves the semidefinite program
    /// \param[out] x The vector minimizing the objective function
    /// \param[in] verbose True to print messages. Default is false
    /// \return The best approximation to the optimal solution
    NT solve(Point& x, bool verbose = false) {
        // initialize
        x = *interiorPoint;
        NT currentMin = objectiveFunction.dot(x);
        int stepsCount = 0;
        NT temperature = diameter;
        NT tempDecreaseFactor = 1.0 - static_cast<NT>(1.0 / std::pow(spectrahedron->dimension(), settings.k));

        // initialize random walk;
        HMC hmc;
        initializeHMC(hmc, diameter);

        // if settings.maxNumSteps is negative there is no
        // bound to the number of steps
        while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) {

            // sample one point with current temperature
            std::list<Point> randPoints;
            hmc.sample(*spectrahedron, x, 1, randPoints);

            // update values;
            x = randPoints.back();
            randPoints.clear();
            currentMin = objectiveFunction.dot(x);
            ++stepsCount;

            if (verbose)
                std::cout << "Step: " << stepsCount << ", Temperature: " << temperature << ", Min: " << currentMin << "\n";

            // decrease the temperature
            temperature *= tempDecreaseFactor;
            hmc.setTemperature(temperature);
        }

        // return the minimum w.r.t. the original objective function
        return currentMin*objectiveFunctionNorm;
    }

    /// Solves the semidefinite program. Stops when reaches a relative error < settings.error
    /// or hits the maximum number of allowed steps
    /// \param[out] x The vector minimizing the objective function
    /// \param[in] verbose True to print messages
    /// \param[in] The exact minimum of the problem
    /// \return The best approximation to the optimal solution
    NT solve(Point& x, bool verbose, const NT exact) {
        // initialize
        x = *interiorPoint;
        NT currentMin = objectiveFunction.dot(x);
        int stepsCount = 0;
        NT temperature = diameter;
        NT tempDecreaseFactor = 1.0 - static_cast<NT>(1.0 / std::pow(spectrahedron->dimension(), settings.k));

        // initialize random walk;
        HMC hmc;
        initializeHMC(hmc, diameter);
        typename HMC::PrecomputedValues hmcPrecomputesValues;

        // if settings.maxNumSteps is negative there is no
        // bound to the number of steps
        while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) {

            // sample one point with current temperature
            std::list<Point> randPoints;

            while (1) {
                hmc.sample(*spectrahedron, x, settings.walkLength, randPoints, hmcPrecomputesValues);

                // if the sampled point is not inside the spectrahedron,
                // get a new one
                if (spectrahedron->isExterior(hmcPrecomputesValues.C)) {
                    if (verbose) std::cout << "Sampled point outside the spectrahedron.\n";
                    randPoints.clear();
                    hmcPrecomputesValues.resetFlags();
                }
                else {
                    // update values;
                    x = randPoints.back();
                    randPoints.clear();
                    break;
                }
            }

            currentMin = objectiveFunction.dot(x);
            ++stepsCount;

            // compute relative error
            NT relError = relativeError(currentMin, exact);

            if (verbose)
                std::cout << "Step: " << stepsCount << ", Temperature: " << temperature << ", Min: " << currentMin
                << ", Relative error: " << relError << "\n";

            // check if we reached desired accuracy
            if (relError < settings.error)
                break;

            // decrease the temperature
            temperature *= tempDecreaseFactor;
            hmc.setTemperature(temperature);
        }

        // return the minimum w.r.t. the original objective function
        return currentMin*objectiveFunctionNorm;
    }

    ~SimulatedAnnealing() {
        // if we computed the interior point
        if (computedInteriorPoint) delete interiorPoint;
    }
};
#endif //VOLESTI_SIMULATEDANNEALING_H
