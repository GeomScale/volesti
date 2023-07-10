// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BOLTZMANN_HMC_WALK_HPP
#define VOLESTI_BOLTZMANN_HMC_WALK_HPP

#include "generators/boost_random_number_generator.hpp"
#include "../sampling/sphere.hpp"

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
struct BoltzmannHMCWalk {
public:

    struct parameters {};
    parameters param;


    /// The implementation of the walk
    /// Currently implemented only for spectrahedra
    /// with template specialization
    ///@tparam ConvexBody a convex body
    ///@tparam RandomNumberGenerator
    template <typename ConvexBody, typename RandomNumberGenerator>
    struct Walk {

        /// The matrix/vector types we use
        typedef typename ConvexBody::PointType Point;
        typedef typename ConvexBody::MT MT;
        typedef typename ConvexBody::VT VT;
        typedef typename Point::FT NT;

        /// A struct containing the parameters for the random walk
        struct Settings {
            /// The number of points to "burn", before keeping the following as a sample
            int walk_length;
            /// For generating random numbers
            RandomNumberGenerator randomNumberGenerator;
            /// The c in the distribution
            VT c;
            /// The T in the distribution
            NT temperature;
            /// The diameter of the body
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
            /// \param[in] diameter The diameter of the convexbody
            /// \param[in] reflectionsBound at each iteration allow reflectionsBound*dimension reflections at most
            /// \param[in] dl approach the boundary with a factor of dl, for numerical stability
            /// \return An instance of this struct
            template <typename Point>
            Settings(const int walkLength, const RandomNumberGenerator &randomNumberGenerator, const Point &c, const NT temperature, const NT diameter,
                     unsigned int reflectionsBound = 100, NT dl = 0.995) : walk_length(walkLength), randomNumberGenerator(randomNumberGenerator),
                                                                          c(c.getCoefficients()),
                                                                          temperature(temperature),
                                                                          diameter(diameter),
                                                                          reflectionsBound(reflectionsBound), dl(dl) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;

        Walk() {}

        /// Constructor
        /// \param[in] settings The settings of the random walk
        Walk(Settings &settings) : settings(settings) {}

        /// Change the settings
        /// \param[in] settings The settings of the random walk
        void setSettings(Settings &settings) {
            this->settings = settings;
        }

        /// Samples random points from the convexbody from the Boltzmann distribution
        /// \param[in] convexbody A convexbody
        /// \param[in] interiorPoint A point in the interior of the convexbody
        /// \param[in] pointsNum The number of points to sample
        /// \param[out] points The list of the sampled points
        /// \tparam Point class Point with NT and VT as declared above in this class
        template <typename Point>
        void apply(ConvexBody &convexbody, Point const & interiorPoint, const unsigned int pointsNum,
                    std::list<Point> &points) {
            // store intermediate results between successive calls of methods
            // of the class convexbody, to avoid repeating computations

            VT p = interiorPoint.getCoefficients();

            // sample #pointsNum points
            for (unsigned int i = 1; i <= pointsNum; ++i) {
                // burn #walk_length points to get one sample
                for (unsigned int j = 0; j < settings.walk_length; ++j) {
                    getNextPoint<Point>(convexbody, p);
                }

                // add the sample in the return list
                points.push_back(Point(p));
            }
        }


        /// A single step of the HMC random walk: choose a direction and walk on the trajectory for a random distance.
        /// If it hits the boundary, the trajectory is reflected. If #reflections < reflectionsBound * dimension, it returns the same point
        /// \param[in] convexbody A convexbody
        /// \param[in, out] p An interior point, and the next point in the random walk
        /// \tparam Point
        template <typename Point>
        void getNextPoint(ConvexBody &convexbody, VT &p) {

            // initialize
            RandomNumberGenerator &rng = settings.randomNumberGenerator;
            boost::random::uniform_real_distribution<> urdist(0, 1);
            const NT dl = settings.dl;
            unsigned int n = convexbody.dimension();
            int reflectionsNum = 0;
            int reflectionsNumBound = settings.reflectionsBound * n;
            VT previousPoint;
            VT p0 = p;

            // choose a distance to walk
            NT T = rng.sample_urdist() * settings.diameter;

            // The trajectory will be of the form a*t^2 + vt + p
            // where a = -c / 2*temperature
            // and at each reflection, v and p will change

            // crate vector a
            VT a = -settings.c / (2 * settings.temperature);

            // The vector v will be a random a direction
            VT v = GetDirection<Point>::apply(n, rng).getCoefficients();

            // Begin a step of the random walk
            // Also, count the reflections and respect the bound
            while (reflectionsNum++ < reflectionsNumBound) {

                // we are at point p and the trajectory a*t^2 + vt + p
                // find how long we can walk till we hit the boundary
                NT lambda = convexbody.positiveQuadIntersection(a, v, p);

                // We just solved a quadratic polynomial eigenvalue system At^2 + Bt + C,
                // where A = lmi(a) - A0, B = lmi(v) - A0 and C = lmi(p)
                // and also did a linearization creating two matrices X, Y.
                // For the subsequent calls, we don't have to compute these matrices from scratch,
                // but can efficiently update them.
                // A remains the same
                // C := A*lambda^2 + B*lambda + C
                // X, Y will be updated in class convexbody
                // Set the flags
                convexbody.set_flags(true);

                // if we can walk the remaining distance without reaching he boundary
                if (T <= lambda) {
                    // set the new point p:= (T^2)*a + T*V + p
                    p += (T * T) * a + T * v;

                    // update matrix C
                    convexbody.update_C(T);
                    return;
                }

                // we hit the boundary and still have to walk
                // don't go all the way to the boundary, for numerical stability
                lambda *= dl;

                // save current and set new point
                previousPoint = p;
                p += (lambda * lambda) * a + lambda * v;

                // update remaining distance we must walk
                T -= lambda;

                // update matrix C
                convexbody.update_C(lambda);
                //precomputedValues.C += (lambda * lambda) * precomputedValues.A + lambda * precomputedValues.B;

                // Set v to have the direction of the trajectory at t = lambda
                // i.e. the gradient of at^2 + vt + p, for t = lambda
                v += (lambda * 2) * a;

                // compute reflected direction
                convexbody.compute_reflection(v, p);
            }

            // if the #reflections exceeded the limit, don't move
            if (reflectionsNum == reflectionsNumBound) {
                p = p0;
                convexbody.set_flags(false);
            }
        }

        /// Sets the temperature in the distribution
        /// \param[in] temperature New value of temperature
        void setTemperature(NT temperature) {
            settings.temperature = temperature;
        }
    };

};

#endif //VOLESTI_BOLTZMANN_HMC_WALK_HPP
