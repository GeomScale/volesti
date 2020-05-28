//
// Created by panagiotis on 2/23/20.
//

#ifndef VOLESTI_HMC_RANDOMWALK_H
#define VOLESTI_HMC_RANDOMWALK_H

#include "spectrahedron.h"
#include "generators/boost_random_number_generator.hpp"
#include "../sampling/sphere.hpp"

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
/// \tparam Point Point type to communicate with the rest library
/// \tparam MT Matrix type for internal computations
/// \tparam VT Vector type for internal computations
/// \tparam RNGType
template<class Point, typename MT, typename VT, typename RNGType>
class HMC_RandomWalk {
public:
    /// The numeric type
    typedef typename Point::FT NT;
    /// The type of the spectrahedron
    typedef Spectrahedron<NT, MT, VT> SPECTRAHEDRON;
    /// Type for internal structure of class Spectrahedron
    typedef typename SPECTRAHEDRON::PrecomputedValues PrecomputedValues;

    /// A struct containing the parameters for the random walk
    struct Settings {
        /// The number of points to "burn", before keeping the following as a sample
        int walk_length;
        /// For generating random numbers
        RNGType rng;
        /// The c in the distribution
        VT c;
        /// The T in the distribution
        NT temperature;
        /// The diameter of the spectrahedron
        NT diameter;
        /// Set the number of allowed reflections at each step: #reflections < reflectionsBound * dimension
        unsigned int reflectionsBound;
        /// When determining we can move d long till we reach the boundary, we walk d*dl, for numerical stability
        NT dl;

        /// Constructs an object of Settings
        /// \param[in] walkLength The number of points to "burn", before keeping the following as a sample
        /// \param[in] rng For generating random numbers
        /// \param[in] c The c in the distribution
        /// \param[in] temperature The T in the distribution
        /// \param[in] diameter The diameter of the spectrahedron
        /// \return An instance of this struct
        Settings(const int walkLength, const RNGType& rng, const Point& c, const NT temperature, const NT diameter,
                 unsigned int reflectionsBound = 10, NT dl = 0.995) : walk_length(walkLength), rng(rng), c(c.getCoefficients()), temperature(temperature),
                 diameter(diameter), reflectionsBound(reflectionsBound), dl(dl) {}

         Settings() {}
    };

    /// The parameters of the random walk
    Settings settings;


    HMC_RandomWalk() {}

    /// Constructor
    /// \param[in] settings The settings of the random walk
    HMC_RandomWalk(Settings& settings) : settings(settings) {}

    /// Change the settings
    /// \param[in] settings The settings of the random walk
    void setSettings(Settings& settings) {
        this->settings = settings;
    }


    /// Samples random points from the spectrahedron from the Boltzmann distribution
    /// \param[in] spectrahedron A spectrahedron
    /// \param[in] interiorPoint A point in the interior of the spectrahedron
    /// \param[in] pointsNum The number of points to sample
    /// \param[out] points The list of the sampled points
    void sample(SPECTRAHEDRON& spectrahedron, const Point& interiorPoint, const unsigned int pointsNum, std::list<Point>& points) {
        // store intermediate results between successive calls of methods
        // of the class spectrahedron, to avoid repeating computations
        PrecomputedValues precomputedValues;
        VT p = interiorPoint.getCoefficients();

        // sample #pointsNum points
        for (unsigned int i = 1; i <= pointsNum; ++i) {
            // burn #walk_length points to get one sample
            for (unsigned int j = 0; j < settings.walk_length; ++j) {
                getNextPoint(spectrahedron, p, precomputedValues);
            }

            // add the sample in the return list
            points.push_back(Point(p));
        }

        // the data in preComputedValues may be out of date in the next call
        precomputedValues.resetFlags();
    }

    /// Samples random points from the spectrahedron from the Boltzmann distribution
    /// \param[in] spectrahedron A spectrahedron
    /// \param[in] interiorPoint A point in the interior of the spectrahedron
    /// \param[in] pointsNum The number of points to sample
    /// \param[out] points The list of the sampled points
    /// \param[in, out] precomputedValues transfer data between sucessive calls
    void sample(SPECTRAHEDRON& spectrahedron, const Point& interiorPoint, const unsigned int pointsNum,
            std::list<Point>& points, PrecomputedValues& precomputedValues) {
        // store intermediate results between successive calls of methods
        // of the class spectrahedron, to avoid repeating computations

        VT p = interiorPoint.getCoefficients();

        // sample #pointsNum points
        for (unsigned int i = 1; i <= pointsNum; ++i) {
            // burn #walk_length points to get one sample
            for (unsigned int j = 0; j < settings.walk_length; ++j) {
                getNextPoint(spectrahedron, p, precomputedValues);
            }

            // add the sample in the return list
            points.push_back(Point(p));
        }

        // the data in preComputedValues may be out of date in the next call
        precomputedValues.resetFlags();
        precomputedValues.computed_C = true;
    }


    // Pick a random direction as a normilized vector
    Point get_direction(const unsigned int dim) {

        boost::normal_distribution<> rdist(0,1);
        NT normal = NT(0);
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);

        Point p(dim);
        NT* data = p.pointerToData();

        //RNGType rng2 = var.rng;
        for (unsigned int i=0; i<dim; i++) {
            *data = rdist(rng);
            normal += *data * *data;
            data++;
        }

        normal=1.0/std::sqrt(normal);
        p *= normal;

        return p;
    }

    /// A single step of the HMC random walk: choose a direction and walk on the trajectory for a random distance.
    /// If it hits the boundary, the trajectory is reflected. If #reflections < reflectionsBound * dimension, it returns the same point
    /// \param[in] spectrahedron A spectrahedron
    /// \param[in, out] p An interior point, and the next point in the random walk
    /// \param[in, out] precomputedValues Data for the methods of the class Spectrahedron
    void getNextPoint(SPECTRAHEDRON& spectrahedron, VT &p, PrecomputedValues& precomputedValues) {

        // initialize
        RNGType &rng = settings.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        const NT dl = settings.dl;
        unsigned int n = spectrahedron.dimension();
        int reflectionsNum = 0;
        int reflectionsNumBound = settings.reflectionsBound * n;
        VT previousPoint;
        VT p0 = p;

        // choose a distance to walk
        NT T = urdist(rng) * settings.diameter;

        // The trajectory will be of the form a*t^2 + vt + p
        // where a = -c / 2*temperature
        // and at each reflection, v and p will change

        // crate vector a
        VT a = - settings.c / (2*settings.temperature);

        // The vector v will be a random a direction
        VT v = get_direction(n).getCoefficients();

        // Begin a step of the random walk
        // Also, count the reflections and respect the bound
        while (reflectionsNum++ < reflectionsNumBound) {

            // we are at point p and the trajectory a*t^2 + vt + p
            // find how long we can walk till we hit the boundary
            NT lambda = spectrahedron.positiveIntersection(a, v, p, precomputedValues);

            // We just solved a quadratic polynomial eigenvalue system At^2 + Bt + C,
            // where A = lmi(a) - A0, B = lmi(v) - A0 and C = lmi(p)
            // and also did a linearization creating two matrices X, Y.
            // For the subsequent calls, we don't have to compute these matrices from scratch,
            // but can efficiently update them.
            // A remains the same
            // C := A*lambda^2 + B*lambda + C
            // X, Y will be updated in class Spectrahedron
            // Set the flags
            precomputedValues.computed_A = true;
            precomputedValues.computed_C = true;
            precomputedValues.computed_XY = true;

            // if we can walk the remaining distance without reaching he boundary
            if (T <= lambda) {
                // set the new point p:= (T^2)*a + T*V + p
                p += (T*T)*a + T*v;

                // update matrix C
                precomputedValues.C += (T*T)*precomputedValues.A + T*precomputedValues.B;
                return;
            }

            // we hit the boundary and still have to walk
            // don't go all the way to the boundary, for numerical stability
            lambda *= dl;

            // save current and set new point
            previousPoint = p;
            p += (lambda*lambda)*a + lambda*v;

            // update remaining distance we must walk
            T -= lambda;

            // update matrix C
            precomputedValues.C += (lambda*lambda)*precomputedValues.A + lambda*precomputedValues.B;

            // Set v to have the direction of the trajectory at t = lambda
            // i.e. the gradient of at^2 + vt + p, for t = lambda
            v += (lambda*2)*a;

            // compute reflected direction
            VT reflectedTrajectory;
            spectrahedron.computeReflection(p, v, reflectedTrajectory, precomputedValues);
            v = reflectedTrajectory;
        }

        // if the #reflections exceeded the limit, don't move
        if(reflectionsNum == reflectionsNumBound)
            p = p0;
    }

    /// Sets the temperature in the distribution
    /// \param[in] temperature New value of temperature
    void setTemperature(NT temperature) {
        settings.temperature = temperature;
    }
};
#endif //VOLESTI_HMC_RANDOMWALK_H
