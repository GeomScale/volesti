// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Angelos Korakitis, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BOLTZMANN_HMC_WALK_HPP
#define VOLESTI_BOLTZMANN_HMC_WALK_HPP

#include "generators/boost_random_number_generator.hpp"
#include "../sampling/sphere.hpp"
#include <limits>
#include <cmath>
#include <algorithm>
#include <vector>

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
/// Includes Robbins-Monro adaptive step sizing for improved performance.

struct BoltzmannHMCWalk {
public:

    struct parameters {};
    parameters param;

    /// The implementation of the walk
    /// Currently implemented only for spectrahedra with template specialization
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
            /// Set the number of allowed reflections at each step
            unsigned int reflectionsBound;
            /// When determining we can move d long till we reach the boundary, we walk d*dl, for numerical stability
            NT dl;
            
            /// Adaptive step size parameters (Robbins-Monro scheme)
            NT target_acceptance_rate;
            NT adaptation_rate;
            int adaptation_window;

            /// Constructs an object of Settings
            template <typename Point>
            Settings(const int walkLength, 
                    const RandomNumberGenerator &randomNumberGenerator, 
                    const Point &c, 
                    const NT temperature, 
                    const NT diameter,
                    unsigned int reflectionsBound = 20, 
                    NT dl = NT(0.995),
                    NT target_acc_rate = NT(0.70),
                    NT adapt_rate = NT(0.2),
                    int adapt_window = 10) 
                : walk_length(walkLength)
                , randomNumberGenerator(randomNumberGenerator)
                , c(c.getCoefficients())
                , temperature(temperature)
                , diameter(diameter)
                , reflectionsBound(reflectionsBound)
                , dl(dl)
                , target_acceptance_rate(target_acc_rate)
                , adaptation_rate(adapt_rate)
                , adaptation_window(adapt_window) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;
        
        /// Robbins-Monro adaptive step sizing state
        NT current_step_scale;
        int total_steps;
        int accepted_steps;
        int adaptation_counter;
        NT last_acceptance_rate;
        
        /// Pre-allocated working vectors for performance
        mutable VT temp_p, temp_v, temp_a;
        mutable VT p_initial_cache, v_initial_cache;

        Walk() : current_step_scale(NT(2.0)), total_steps(0), accepted_steps(0), 
                 adaptation_counter(0), last_acceptance_rate(NT(0.5)) {}

        /// Constructor
        Walk(Settings &settings) : settings(settings), current_step_scale(NT(2.0)), 
                                  total_steps(0), accepted_steps(0), adaptation_counter(0),
                                  last_acceptance_rate(NT(0.5)) {}

        /// Change the settings
        void setSettings(Settings &settings) {
            this->settings = settings;
        }

        /// Samples random points from the convexbody from the Boltzmann distribution
        template <typename Point>
        void apply(ConvexBody &convexbody, Point const & interiorPoint, const unsigned int pointsNum,
                    std::list<Point> &points) {
            
            VT p = interiorPoint.getCoefficients();
            convexbody.resetFlags();

            // Sample #pointsNum points
            for (unsigned int i = 1; i <= pointsNum; ++i) {
                // Burn #walk_length points to get one sample
                for (unsigned int j = 0; j < settings.walk_length; ++j) {
                    getNextPoint<Point>(convexbody, p);
                }
                // Add the sample in the return list
                points.push_back(Point(p));
            }
        }

        /// A single step of the HMC random walk with Metropolis acceptance and adaptive step sizing
        template <typename Point>
        void getNextPoint(ConvexBody &convexbody, VT &p) {

            RandomNumberGenerator &rng = settings.randomNumberGenerator;
            const NT dl = settings.dl;
            const unsigned int n = convexbody.dimension();
            const int reflectionsNumBound = settings.reflectionsBound * n;
            
            // Store initial state - reuse cache vectors
            p_initial_cache = p;
            temp_p = p; // Current position during trajectory
            
            // Adaptive trajectory length
            NT T = rng.sample_urdist() * settings.diameter * current_step_scale;
            
            // Early termination for very small steps
            if (T < NT(1e-8)) {
                recordAcceptance(false);
                return;
            }

            // Pre-compute acceleration: a = -c / (2*temperature)
            const NT temp_inv = NT(1.0) / (NT(2.0) * settings.temperature);
            temp_a = settings.c;
            temp_a *= -temp_inv; // In-place multiply

            // Sample initial momentum - reuse cache
            v_initial_cache = GetDirection<Point>::apply(n, rng).getCoefficients();
            temp_v = v_initial_cache; // Current velocity

            // Pre-compute initial energy components
            const NT U_initial = settings.c.dot(p_initial_cache) / settings.temperature;
            const NT K_initial = NT(0.5) * v_initial_cache.squaredNorm();
            const NT H_initial = U_initial + K_initial;

            convexbody.resetFlags();
            int reflectionsNum = 0;

            // Main trajectory loop
            while (reflectionsNum < reflectionsNumBound && T > NT(1e-8)) {

                NT lambda;
                try {
                    lambda = convexbody.positiveQuadIntersection(temp_a, temp_v, temp_p);
                } catch (...) {
                    break;
                }

                if (lambda <= NT(0) || !std::isfinite(lambda)) break;

                convexbody.set_flags(true);

                // Complete trajectory without boundary hit
                if (T <= lambda) {
                    // In-place position update: p += T*v + T^2*a
                    temp_p += T * temp_v;
                    temp_p += (T * T) * temp_a;
                    
                    // In-place velocity update: v += T*a
                    temp_v += T * temp_a;
                    
                    convexbody.update_C(T);
                    
                    // Movement check - reuse computation
                    temp_p -= p_initial_cache; // temp_p now contains displacement
                    NT movement_sq = temp_p.squaredNorm();
                    temp_p += p_initial_cache; // restore temp_p
                    
                    if (movement_sq < NT(1e-12)) {
                        recordAcceptance(false);
                        convexbody.set_flags(false);
                        return;
                    }
                    
                    // Final energy computation
                    NT U_final = settings.c.dot(temp_p) / settings.temperature;
                    NT K_final = NT(0.5) * temp_v.squaredNorm();
                    NT delta_H = (U_final + K_final) - H_initial;
                    
                    // Metropolis test
                    bool accept = (delta_H <= NT(0)) || (rng.sample_urdist() < std::exp(-delta_H));
                    
                    if (accept) {
                        p = temp_p;
                        recordAcceptance(true);
                    } else {
                        recordAcceptance(false);
                        convexbody.set_flags(false);
                    }
                    return;
                }

                // Hit boundary - in-place updates
                lambda *= dl;
                
                // Update position: p += λ*v + λ²*a
                temp_p += lambda * temp_v;
                temp_p += (lambda * lambda) * temp_a;
                
                T -= lambda;
                convexbody.update_C(lambda);
                
                // Update velocity: v += 2λ*a
                temp_v += (lambda * NT(2.0)) * temp_a;

                // Reflection
                Point v_point(temp_v);
                Point p_point(temp_p);

                try {
                    convexbody.compute_reflection(v_point, p_point);
                    temp_v = v_point.getCoefficients();
                } catch (...) {
                    break;
                }

                reflectionsNum++;
            }

            // Trajectory failed
            recordAcceptance(false);
            convexbody.set_flags(false);
        }

        /// Sets the temperature in the distribution
        void setTemperature(NT temperature) {
            settings.temperature = temperature;
        }
        
        /// Get current acceptance rate for monitoring
        NT getAcceptanceRate() const {
            if (total_steps == 0) return last_acceptance_rate;
            return static_cast<NT>(accepted_steps) / total_steps;
        }

    private:
        
        /// Record acceptance and apply simple adaptive step sizing
        void recordAcceptance(bool accepted) {
            if (accepted) accepted_steps++;
            total_steps++;
            adaptation_counter++;
            
            // Apply adaptation every adaptation_window steps
            if (adaptation_counter >= settings.adaptation_window) {
                NT acc_rate = static_cast<NT>(accepted_steps) / total_steps;
                last_acceptance_rate = acc_rate;
                NT target = settings.target_acceptance_rate;
                
                // Aggressive adaptation for extreme acceptance rates
                if (acc_rate > NT(0.95)) {
                    current_step_scale *= NT(1.5);  // Very high -> big increase
                } else if (acc_rate > target + NT(0.15)) {
                    current_step_scale *= NT(1.2);  // High -> increase
                } else if (acc_rate < NT(0.2)) {
                    current_step_scale *= NT(0.4);  // Very low -> big decrease
                } else if (acc_rate < target - NT(0.15)) {
                    current_step_scale *= NT(0.7);  // Low -> decrease
                }
                
                // Clamp step scale to reasonable bounds
                current_step_scale = std::max(NT(0.01), std::min(current_step_scale, NT(50.0)));
                
                // Reset all counters for next adaptation window
                adaptation_counter = 0;
                accepted_steps = 0;
                total_steps = 0;
            }
        }
    };
};

#endif //VOLESTI_BOLTZMANN_HMC_WALK_HPP