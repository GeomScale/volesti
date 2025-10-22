// VolEsti (volume computation and sampling library)

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BOLTZMANN_HMC_WALK_HPP
#define VOLESTI_BOLTZMANN_HMC_WALK_HPP

#include "generators/boost_random_number_generator.hpp"
#include "../sampling/sphere.hpp"
#include <algorithm>
#include <cmath>
#include <list>
#include <stdexcept>
#include <vector>

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
/// This implementation includes adaptive step sizing, velocity persistence, and temperature scaling
/// for improved sampling efficiency, particularly useful for simulated annealing.
struct BoltzmannHMCWalk {
public:
    struct parameters {};
    parameters param;

    /// The implementation of the walk
    /// \tparam ConvexBody A convex body
    /// \tparam RandomNumberGenerator A random number generator
    template <typename ConvexBody, typename RandomNumberGenerator>
    struct Walk {
        /// The matrix/vector types we use
        typedef typename ConvexBody::PointType Point;
        typedef typename ConvexBody::MT MT;
        typedef typename ConvexBody::VT VT;
        typedef typename Point::FT NT;

        /// A struct containing the parameters for the random walk
        struct Settings {
            /// The number of points to "burn" before keeping the following as a sample
            int walk_length;
            
            /// For generating random numbers
            RandomNumberGenerator randomNumberGenerator;
            
            /// The c in the Boltzmann distribution e^(-c*x/T)
            VT c;
            
            /// The temperature T in the distribution
            NT temperature;
            
            /// The diameter of the convex body
            NT diameter;
            
            /// Set the maximum number of allowed reflections at each step
            unsigned int reflectionsBound;
            
            /// Initial step size scale factor (will be adapted during sampling)
            NT initial_step_scale;

            /// Constructs an object of Settings
            /// \param[in] walkLength The number of points to "burn" before keeping the following as a sample
            /// \param[in] randomNumberGenerator For generating random numbers
            /// \param[in] c The c in the distribution e^(-c*x/T)
            /// \param[in] temperature The T in the distribution
            /// \param[in] diameter The diameter of the convex body
            /// \param[in] reflectionsBound Maximum number of reflections allowed per step
            /// \param[in] stepScale Initial step size scaling factor (default 1.0)
            /// \return An instance of this struct
            template <typename Point>
            Settings(const int walkLength,
                    const RandomNumberGenerator &randomNumberGenerator,
                    const Point &c,
                    const NT temperature,
                    const NT diameter,
                    unsigned int reflectionsBound = 2000,
                    NT stepScale = NT(1.0))
                : walk_length(walkLength)
                , randomNumberGenerator(randomNumberGenerator)
                , c(c.getCoefficients())
                , temperature(temperature)
                , diameter(diameter)
                , reflectionsBound(reflectionsBound)
                , initial_step_scale(stepScale) {
                // Validate parameters
                if (temperature <= NT(0)) {
                    throw std::invalid_argument("Temperature must be positive");
                }
                if (diameter <= NT(0)) {
                    throw std::invalid_argument("Diameter must be positive");
                }
            }

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;
        
        /// Adaptive step size control: current step scale factor
        NT current_step_scale;
        
        /// Reference temperature for relative scaling (used in simulated annealing)
        NT reference_temperature;
        
        /// Last temperature value to detect temperature changes
        NT last_temperature;
        
        /// Statistics for adaptive step sizing
        unsigned int total_steps;
        unsigned int accepted_steps;
        
        /// Flag to enable/disable adaptive step sizing
        bool adapt_enabled;

        /// Working buffers to avoid repeated allocations
        mutable VT position_buffer;
        mutable VT velocity_buffer;
        mutable VT acceleration_buffer;
        
        /// Previous velocity for velocity persistence (improves mixing)
        mutable VT previous_velocity;
        mutable bool has_previous_velocity;

        /// Default constructor
        Walk() : current_step_scale(NT(1.0)), 
                 reference_temperature(NT(1.0)),
                 last_temperature(NT(-1)),
                 total_steps(0), 
                 accepted_steps(0), 
                 adapt_enabled(true),
                 has_previous_velocity(false) {}

        /// Constructor with settings
        /// \param[in] settings The settings of the random walk
        Walk(Settings &settings) : settings(settings),
                                   current_step_scale(settings.initial_step_scale),
                                   reference_temperature(settings.temperature),
                                   last_temperature(NT(-1)),
                                   total_steps(0),
                                   accepted_steps(0),
                                   adapt_enabled(true),
                                   has_previous_velocity(false) {}

        /// Change the settings
        /// \param[in] settings The new settings of the random walk
        void setSettings(Settings &settings) {
            this->settings = settings;
        }

        /// Samples random points from the convex body using the Boltzmann distribution
        /// \param[in] convexbody A convex body
        /// \param[in] interiorPoint A point in the interior of the convex body
        /// \param[in] pointsNum The number of points to sample
        /// \param[out] points The list of the sampled points
        /// \tparam Point class Point with NT and VT as declared above in this class
        template <typename Point>
        void apply(ConvexBody &convexbody, Point const & interiorPoint, 
                   const unsigned int pointsNum, std::list<Point> &points) {
            
            // Get the dimension and initialize working buffers
            const unsigned int dim = convexbody.dimension();
            if (position_buffer.size() != dim) {
                position_buffer.resize(dim);
                velocity_buffer.resize(dim);
                acceleration_buffer.resize(dim);
                previous_velocity.resize(dim);
            }

            // Start from the interior point
            VT p = interiorPoint.getCoefficients();

            // Sample #pointsNum points
            for (unsigned int i = 1; i <= pointsNum; ++i) {
                // Burn #walk_length points to get one sample
                for (unsigned int j = 0; j < settings.walk_length; ++j) {
                    getNextPoint<Point>(convexbody, p);
                }
                // Add the sample to the return list
                points.push_back(Point(p));
            }
        }

        /// A single step of the HMC random walk with Metropolis-Hastings acceptance.
        /// The trajectory follows the Hamiltonian dynamics H = U(x) + K(v), where
        /// U(x) = c*x/T (potential energy) and K(v) = 0.5*||v||^2 (kinetic energy).
        /// 
        /// \param[in] convexbody A convex body
        /// \param[in,out] p An interior point, updated to the next point in the random walk
        /// \tparam Point
        template <typename Point>
        void getNextPoint(ConvexBody &convexbody, VT &p) {
            // Reset body flags when temperature changes (important for simulated annealing)
            if (last_temperature != settings.temperature) {
                resetBodyFlags(convexbody);
                last_temperature = settings.temperature;
            }

            // Start with current position as proposal
            VT proposed_position = p;
            
            // Compute trajectory length with temperature scaling
            // Scale relative to reference temperature: higher T → longer steps
            NT temp_scale = std::sqrt(settings.temperature / reference_temperature);
            temp_scale = std::max(temp_scale, NT(0.01)); // Minimum scale to avoid numerical issues
            NT trajectory_length = settings.diameter * current_step_scale * temp_scale;
            
            // Compute acceleration: a = -c/T (gradient of potential energy)
            // In Hamiltonian dynamics, acceleration comes from the gradient of U(x)
            const NT invT = NT(1) / settings.temperature;
            acceleration_buffer = (-invT) * settings.c;

            // Sample initial velocity with persistence for better mixing
            // Velocity persistence: v_new = rho * v_old + sqrt(1-rho^2) * noise
            // This creates correlation between steps, improving exploration
            RandomNumberGenerator &rng = settings.randomNumberGenerator;
            const NT rho = NT(0.9);  // Persistence parameter (0 = no persistence, 1 = full persistence)
            const NT noise = std::sqrt(NT(1) - rho * rho);
            const unsigned int dim = position_buffer.size();
            
            if (has_previous_velocity && previous_velocity.size() == dim) {
                // Use persistent velocity
                for (unsigned i = 0; i < dim; ++i) {
                    velocity_buffer[i] = rho * previous_velocity[i] + 
                                        noise * rng.sample_ndist();
                }
            } else {
                // First step or dimension changed: sample fresh velocity
                for (unsigned i = 0; i < dim; ++i) {
                    velocity_buffer[i] = rng.sample_ndist();
                }
            }

            // Store initial velocity and compute initial Hamiltonian
            const VT v0 = velocity_buffer;
            const NT H0 = proposed_position.dot(settings.c) * invT +  // Potential energy U(x)
                         NT(0.5) * v0.squaredNorm();                   // Kinetic energy K(v)

            // Execute the Hamiltonian trajectory
            unsigned int reflections = 0;
            const NT used = executeTrajectory(convexbody, proposed_position, 
                                             trajectory_length, reflections);

            // Compute final Hamiltonian for Metropolis-Hastings acceptance test
            const NT H1 = proposed_position.dot(settings.c) * invT + 
                         NT(0.5) * velocity_buffer.squaredNorm();
            const NT dH = H1 - H0;  // Change in Hamiltonian (should be near zero ideally)

            // Metropolis-Hastings acceptance criterion: accept with probability min(1, exp(-dH))
            bool accepted = false;
            if (std::isfinite(dH)) {
                if (dH <= NT(0)) {
                    // Hamiltonian decreased: always accept
                    accepted = true;
                } else if (dH > NT(50)) {
                    // exp(-50) ≈ 2e-22, effectively zero probability
                    accepted = false;
                } else {
                    // Probabilistic acceptance based on exp(-dH)
                    const NT u = std::max(NT(1e-300), 
                                         static_cast<NT>(rng.sample_urdist()));
                    accepted = (-std::log(u)) > dH;
                }
            }

            // Update position and velocity based on acceptance
            if (accepted) {
                p = proposed_position;
                previous_velocity = velocity_buffer;
                has_previous_velocity = true;
                accepted_steps++;
            } else {
                // Reject: keep current position, but update velocity (for persistence)
                velocity_buffer = v0;
                previous_velocity = velocity_buffer;
                has_previous_velocity = true;
            }

            total_steps++;
            
            // Adaptive step size adjustment 
            // Goal: maximize trajectory completion while minimizing reflections
            if (adapt_enabled) {
                // Measure trajectory completion ratio
                const NT completion = (trajectory_length > NT(1e-12)) ? 
                                     (used / trajectory_length) : NT(0);
                
                // Measure reflection rate (normalized by distance traveled)
                const NT norm_dist = std::max(NT(1), used / settings.diameter);
                const NT reflection_rate = NT(reflections) / norm_dist;

                // Define target ranges
                const NT good_completion = NT(0.85);    // High completion is good
                const NT poor_completion = NT(0.45);    // Low completion is bad
                const NT high_reflections = NT(0.7);    // Too many reflections is bad - those slow down the walk

                // Compute step size adjustment
                NT adjustment = NT(0);
                if (completion > good_completion && reflection_rate < high_reflections) {
                    // Good trajectory: increase step size
                    adjustment = NT(0.05);
                } else if (completion < poor_completion || reflection_rate > high_reflections) {
                    // Poor trajectory: decrease step size
                    adjustment = -NT(0.05);
                } else {
                    // Moderate trajectory: adjust toward target
                    const NT target = (good_completion + poor_completion) / NT(2);
                    adjustment = NT(0.02) * (completion - target);
                }

                // Apply adjustment (exponential update to ensure positivity)
                current_step_scale *= std::exp(adjustment);
                // Clamp to reasonable range
                current_step_scale = std::max(NT(0.1), std::min(current_step_scale, NT(1e8)));
            }
        }

        /// Sets the temperature in the distribution
        /// \param[in] temperature New value of temperature (must be positive)
        void setTemperature(NT temperature) {
            if (temperature <= NT(0)) {
                throw std::invalid_argument("Temperature must be positive");
            }
            settings.temperature = temperature;
        }

        /// Gets the current effective step size (including temperature scaling)
        /// \return The current step size
        NT getStepSize() const {
            NT temp_scale = std::sqrt(settings.temperature / reference_temperature);
            temp_scale = std::max(temp_scale, NT(0.1));
            return settings.diameter * current_step_scale * temp_scale;
        }

        /// Gets the current acceptance rate of the Metropolis-Hastings test
        /// \return The acceptance rate (between 0 and 1)
        NT getAcceptanceRate() const {
            return (total_steps == 0) ? NT(0.5) : 
                   static_cast<NT>(accepted_steps) / static_cast<NT>(total_steps);
        }

        /// Enable or disable adaptive step sizing
        /// \param[in] enabled True to enable, false to disable
        void setAdaptationEnabled(bool enabled) { 
            adapt_enabled = enabled; 
        }

        /// Manually set the step scale factor
        /// \param[in] scale The new step scale factor
        void setStepScale(NT scale) { 
            current_step_scale = scale; 
        }

        /// Set the maximum number of reflections allowed per trajectory
        /// \param[in] maxReflections The maximum number of reflections
        void setMaxReflections(unsigned int maxReflections) {
            settings.reflectionsBound = maxReflections;
        }

    private:
        /// Reset body computation cache on temperature change (SFINAE pattern)
        /// This is important when the body caches computations that depend on temperature
        template <typename Body>
        static auto resetBodyFlags(Body& b) -> decltype(b.resetFlags(), void()) {
            b.resetFlags();
        }
        
        /// Fallback for bodies that don't have resetFlags method
        template <typename Body>
        static void resetBodyFlags(...) {}

        /// Helper function: compute position along parabolic trajectory
        /// Trajectory is x(t) = x0 + v*t + 0.5*a*t^2
        /// \param[in] pos Initial position x0
        /// \param[in] vel Initial velocity v
        /// \param[in] acc Acceleration a
        /// \param[in] time Time parameter t
        /// \return Position at time t
        inline VT computeTrajectoryPoint(const VT& pos, const VT& vel, 
                                         const VT& acc, NT time) const {
            return pos + time * vel + (NT(0.5) * time * time) * acc;
        }

        /// Execute a Hamiltonian trajectory with boundary reflections and epsilon-inward handling.
        /// The trajectory follows x(t) = x0 + v*t + 0.5*a*t^2 until hitting the boundary,
        /// at which point the velocity is reflected.
        /// 
        /// \param[in] body The convex body
        /// \param[in,out] position Starting position, updated to final position
        /// \param[in] max_length Maximum trajectory length to execute
        /// \param[out] reflections Number of reflections that occurred
        /// \return Actual distance traveled along trajectory
        NT executeTrajectory(ConvexBody& body, VT& position,
                            NT max_length, unsigned int& reflections) {
            reflections = 0;
            NT remaining = max_length;
            const NT eps = NT(1e-10);  // Safety margin from boundary
            const unsigned int maxRefl = settings.reflectionsBound;

            // Continue until trajectory is complete or we exceed reflection limit
            while (remaining > NT(1e-12)) {
                // Early exit if too many reflections occurred
                if (reflections >= maxRefl) break;

                // Find how long we can walk along the quadratic trajectory
                // before hitting the boundary (solves a quadratic eigenvalue problem)
                NT t = NT(-1);
                try {
                    t = body.positiveQuadIntersection(acceleration_buffer, 
                                                     velocity_buffer, position);
                } catch (...) {
                    // Numerical issue: stop trajectory
                    break;
                }

                // Case 1: Can complete the trajectory without hitting boundary
                if (t >= remaining - NT(1e-12)) {
                    VT trial = computeTrajectoryPoint(position, velocity_buffer, 
                                                      acceleration_buffer, remaining);

                    // Safety check: verify we're still inside
                    // (numerical errors can sometimes push us outside)
                    if (body.isExterior(trial)) {
                        // Try going halfway instead
                        const NT half = remaining * NT(0.5);
                        trial = computeTrajectoryPoint(position, velocity_buffer, 
                                                      acceleration_buffer, half);
                        if (body.isExterior(trial)) break;  // Give up if still outside
                        remaining = half;
                    }

                    // Update position and velocity
                    position = trial;
                    velocity_buffer += remaining * acceleration_buffer;
                    body.update_C(remaining);
                    return max_length;  // Successfully completed full trajectory
                }

                // Case 2: Will hit boundary before completing trajectory
                // Move to epsilon-inward from boundary for numerical stability
                NT safe_t = (t > eps) ? (t - eps) : (t * NT(0.5));
                
                VT boundary_pos = computeTrajectoryPoint(position, velocity_buffer, 
                                                         acceleration_buffer, safe_t);

                // Fallback: if epsilon-inward point is outside, try half distance
                if (body.isExterior(boundary_pos)) {
                    safe_t = t * NT(0.5);
                    boundary_pos = computeTrajectoryPoint(position, velocity_buffer, 
                                                         acceleration_buffer, safe_t);
                    if (body.isExterior(boundary_pos)) break;  // Give up if still outside
                }

                // Update position and velocity to boundary point
                position = boundary_pos;
                velocity_buffer += safe_t * acceleration_buffer;
                body.update_C(safe_t);
                remaining -= safe_t;

                // Reflect velocity at the boundary
                // The body computes the normal and reflects the velocity
                Point v_pt(velocity_buffer);
                Point x_pt(position);
                try {
                    body.compute_reflection(v_pt, x_pt);
                    velocity_buffer = v_pt.getCoefficients();
                    reflections++;
                } catch (...) {
                    // Reflection failed: stop trajectory
                    break;
                }
            }

            // Return actual distance traveled (may be less than max_length)
            return max_length - remaining;
        }
    };
};

#endif // VOLESTI_BOLTZMANN_HMC_WALK_HPP