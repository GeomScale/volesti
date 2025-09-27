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
#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
/// Implementation with temperature-aware adaptive step sizing for convex optimization.

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
            /// The c in the distribution (linear coefficients)
            VT c;
            /// The T in the distribution (temperature)
            NT temperature;
            /// The diameter of the body
            NT diameter;
            /// Set the number of allowed reflections at each step
            unsigned int reflectionsBound;
            /// Initial step scale factor
            NT initial_step_scale;

            /// Constructs an object of Settings
            template <typename Point>
            Settings(const int walkLength,
                    const RandomNumberGenerator &randomNumberGenerator,
                    const Point &c,
                    const NT temperature,
                    const NT diameter,
                    unsigned int reflectionsBound = 200,
                    NT stepScale = NT(1.0))
                : walk_length(walkLength)
                , randomNumberGenerator(randomNumberGenerator)
                , c(c.getCoefficients())
                , temperature(temperature)
                , diameter(diameter)
                , reflectionsBound(reflectionsBound)
                , initial_step_scale(stepScale) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;

        /// Adaptive step sizing state
        NT current_step_scale;
        NT reference_temperature;
        NT last_temperature;
        unsigned int total_steps;
        unsigned int accepted_steps;
        bool adapt_enabled;

        /// Pre-allocated working vectors for performance
        mutable VT position_buffer, velocity_buffer, acceleration_buffer;
        mutable VT previous_velocity;
        mutable bool has_previous_velocity;

        /// Constructor
        Walk() : current_step_scale(NT(1.0)), reference_temperature(NT(1.0)),
                 last_temperature(NT(-1)), total_steps(0), accepted_steps(0),
                 adapt_enabled(true), has_previous_velocity(false) {}

        /// Constructor with settings
        Walk(Settings &settings) : settings(settings),
                                   current_step_scale(settings.initial_step_scale),
                                   reference_temperature(settings.temperature),
                                   last_temperature(NT(-1)),
                                   total_steps(0),
                                   accepted_steps(0),
                                   adapt_enabled(true),
                                   has_previous_velocity(false) {}

        /// Change the settings
        void setSettings(Settings &settings) {
            this->settings = settings;
        }

        /// Samples random points from the convexbody from the Boltzmann distribution
        template <typename Point>
        void apply(ConvexBody &convexbody, Point const & interiorPoint, const unsigned int pointsNum,
                   std::list<Point> &points) {

            validateConfiguration();
            initializeBuffers(convexbody.dimension());
            initializeStepScale(convexbody.dimension());

            VT p = interiorPoint.getCoefficients();

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

            maybeResetForTemperature(convexbody);

            VT proposed_position = p;
            const NT trajectory_length = computeTrajectoryLength();
            if (trajectory_length <= NT(1e-12)) return;

            computeAcceleration();
            initializeVelocity(convexbody.dimension());

            const VT v0 = velocity_buffer;
            const NT H0 = computeHamiltonian(proposed_position, v0);

            unsigned int reflections = 0;
            const NT used = executeTrajectory(convexbody, proposed_position, trajectory_length, reflections);

            const bool accepted = metropolisAccept(H0, proposed_position);
            if (accepted) {
                p = proposed_position;
                previous_velocity = velocity_buffer;
                has_previous_velocity = true;
                accepted_steps++;
            } else {
                velocity_buffer = v0;
                previous_velocity = velocity_buffer;
                has_previous_velocity = true;
            }

            total_steps++;
            adaptStepSize(trajectory_length, used, reflections);
        }

        /// Sets the temperature in the distribution
        void setTemperature(NT temperature) {
            if (temperature <= NT(0)) throw std::invalid_argument("Temperature must be positive");
            settings.temperature = temperature;
        }

        /// Get current step size
        NT getStepSize() const { return computeTrajectoryLength(); }

        /// Get current acceptance rate for monitoring
        NT getAcceptanceRate() const {
            if (total_steps == 0) return NT(0.5);
            return static_cast<NT>(accepted_steps) / static_cast<NT>(total_steps);
        }

        /// Enable or disable adaptation
        void setAdaptationEnabled(bool enabled) { adapt_enabled = enabled; }

        /// Set step scale directly
        void setStepScale(NT scale) { current_step_scale = scale; }

        /// Set maximum reflections
        void setMaxReflections(unsigned int maxReflections) {
            settings.reflectionsBound = maxReflections;
        }

    private:

        /// Minimum and maximum scale factors
        static constexpr NT kMinScale() { return NT(0.1); }
        static constexpr NT kMaxScale() { return NT(1e4); }
        static constexpr NT kEps() { return NT(1e-12); }

        /// Clamp a value between bounds
        template <typename T>
        static T clamp(const T& x, const T& lo, const T& hi) {
            return (x < lo) ? lo : ((x > hi) ? hi : x);
        }

        /// Cache reset when temperature changes
        template <typename Body>
        static auto resetFlagsIfSupported(Body& b, int) -> decltype(b.resetFlags(), void()) {
            b.resetFlags();
        }
        static void resetFlagsIfSupported(...) {}

        template <typename Body>
        void maybeResetForTemperature(Body& body) {
            if (last_temperature != settings.temperature) {
                resetFlagsIfSupported(body, 0);
                last_temperature = settings.temperature;
            }
        }

        /// Validate configuration parameters
        void validateConfiguration() const {
            if (settings.temperature <= NT(0)) {
                throw std::invalid_argument("Temperature must be positive");
            }
            if (settings.diameter <= NT(0)) {
                throw std::invalid_argument("Body diameter must be positive");
            }
        }

        /// Initialize working buffers
        void initializeBuffers(unsigned int dimension) {
            if (position_buffer.size() != dimension) {
                position_buffer.resize(dimension);
                velocity_buffer.resize(dimension);
                acceleration_buffer.resize(dimension);
            }
        }

        /// Initialize step scale for high dimensions
        void initializeStepScale(unsigned int dimension) {
            // For high-dimensional convex optimization: start with larger steps
            if (dimension >= 100 && current_step_scale < NT(2)) {
                current_step_scale = std::sqrt(NT(dimension) / NT(25));
            }

            // Increase reflection budget for high dimensions
            if (dimension >= 100 && settings.reflectionsBound < 300) {
                settings.reflectionsBound = std::min(unsigned(1000),
                    unsigned(10 * std::sqrt(static_cast<NT>(dimension))));
            }
        }

        /// Compute trajectory length: CRITICAL for low temperature
        /// When T approaches 0, acceleration magnitude |a| = |c|/T approaches infinity
        /// Must scale steps down with sqrt(T) to maintain reasonable trajectory curvature
        NT computeTrajectoryLength() const {
            NT temp_scale = std::sqrt(settings.temperature / reference_temperature);
            temp_scale = std::max(temp_scale, NT(0.01));  // Never below 1% of reference

            NT length = settings.diameter * current_step_scale * temp_scale;
            NT min_length = settings.diameter * NT(1e-8);
            return std::max(length, min_length);
        }

        /// Compute acceleration from linear potential
        void computeAcceleration() {
            const NT invT = NT(1) / settings.temperature;
            acceleration_buffer = (-invT) * settings.c;
        }

        /// Initialize velocity with high persistence for convex optimization
        void initializeVelocity(unsigned int dimension) {
            RandomNumberGenerator &rng = settings.randomNumberGenerator;
            // High persistence for convex optimization (not exploration)
            const NT rho = (dimension >= 100) ? NT(0.95) : NT(0.90);
            const NT noise = std::sqrt(NT(1) - rho * rho);

            if (has_previous_velocity && previous_velocity.size() == dimension) {
                for (unsigned i = 0; i < dimension; ++i) {
                    velocity_buffer[i] = rho * previous_velocity[i] +
                                        noise * rng.sample_ndist();
                }
            } else {
                for (unsigned i = 0; i < dimension; ++i) {
                    velocity_buffer[i] = rng.sample_ndist();
                }
            }
        }

        /// Compute Hamiltonian (total energy)
        NT computeHamiltonian(const VT& x, const VT& v) const {
            const NT invT = NT(1) / settings.temperature;
            return x.dot(settings.c) * invT + NT(0.5) * v.squaredNorm();
        }

        /// Outcome of a trajectory step
        struct StepOutcome {
            bool succeeded;
            bool hit_boundary;
            NT distance_covered;
        };

        /// Execute trajectory with reflections
        NT executeTrajectory(ConvexBody& body, VT& position,
                            NT max_length, unsigned int& reflections) {
            reflections = 0;
            NT remaining = max_length;
            unsigned int consecutive_grazes = 0;

            while (reflections < settings.reflectionsBound && remaining > kEps()) {
                // Find intersection time with boundary
                NT t = findIntersectionTime(body, position);

                if (!(t > NT(0)) || !std::isfinite(static_cast<double>(t))) {
                    // Fallback: conservative advance
                    StepOutcome s = conservativeAdvance(body, position, remaining);
                    if (!s.succeeded) break;
                    remaining -= s.distance_covered;
                    continue;
                }

                if (t >= remaining - kEps()) {
                    // Can complete trajectory without hitting boundary
                    StepOutcome s = completeTrajectory(body, position, remaining);
                    if (!s.succeeded) break;
                    remaining -= s.distance_covered;
                    continue;
                }

                // Hit boundary before trajectory end
                StepOutcome s = advanceToBoundary(body, position, t);
                if (!s.succeeded) break;
                remaining -= s.distance_covered;

                // Detect grazing (hitting boundary with near-zero progress)
                if (s.distance_covered <= max_length * NT(1e-10)) {
                    if (++consecutive_grazes >= 3) break;
                } else {
                    consecutive_grazes = 0;
                }

                // Reflect velocity at boundary
                if (!reflectAtBoundary(body, position)) break;
                reflections++;

                // Nudge slightly inside after reflection
                if (!nudgeInside(body, position, remaining * NT(1e-8))) break;
            }

            return max_length - remaining;
        }

        /// Find time to boundary intersection
        NT findIntersectionTime(ConvexBody& body, const VT& position) const {
            try {
                return body.positiveQuadIntersection(acceleration_buffer, velocity_buffer, position);
            } catch (...) {
                return std::numeric_limits<NT>::quiet_NaN();
            }
        }

        /// Conservative advance when intersection fails
        StepOutcome conservativeAdvance(ConvexBody& body, VT& position, NT max_dist) {
            NT step = std::min(max_dist, computeTrajectoryLength() * NT(0.1));
            const NT min_step = settings.diameter * NT(1e-12);

            for (int it = 0; it < 40 && step > min_step; ++it) {
                VT trial = position + step * velocity_buffer +
                          NT(0.5) * step * step * acceleration_buffer;
                if (!body.isExterior(trial)) {
                    position = trial;
                    velocity_buffer += step * acceleration_buffer;
                    body.update_C(step);
                    return {true, false, step};
                }
                step *= NT(0.5);
            }
            return {false, false, NT(0)};
        }

        /// Complete trajectory without boundary hit
        StepOutcome completeTrajectory(ConvexBody& body, VT& position, NT dist) {
            VT trial = position + dist * velocity_buffer +
                      NT(0.5) * dist * dist * acceleration_buffer;

            if (body.isExterior(trial)) {
                // Try half-step as fallback
                NT half = dist * NT(0.5);
                trial = position + half * velocity_buffer +
                       NT(0.5) * half * half * acceleration_buffer;
                if (body.isExterior(trial)) return {false, false, half};
                dist = half;
            }

            position = trial;
            velocity_buffer += dist * acceleration_buffer;
            body.update_C(dist);
            return {true, false, dist};
        }

        /// Advance to boundary
        StepOutcome advanceToBoundary(ConvexBody& body, VT& position, NT t) {
            VT bnd = position + t * velocity_buffer +
                    NT(0.5) * t * t * acceleration_buffer;

            if (body.isExterior(bnd)) {
                // Numerical overshoot: try half-step
                NT half = t * NT(0.5);
                bnd = position + half * velocity_buffer +
                     NT(0.5) * half * half * acceleration_buffer;
                if (body.isExterior(bnd)) return {false, true, half};
                t = half;
            }

            position = bnd;
            velocity_buffer += t * acceleration_buffer;
            body.update_C(t);
            return {true, true, t};
        }

        /// Reflect velocity at boundary
        bool reflectAtBoundary(ConvexBody& body, const VT& pos_at_boundary) {
            Point v_pt(velocity_buffer);
            Point x_pt(pos_at_boundary);
            try {
                body.compute_reflection(v_pt, x_pt);
                velocity_buffer = v_pt.getCoefficients();
                return true;
            } catch (...) {
                return false;
            }
        }

        /// Nudge position slightly inside after reflection
        bool nudgeInside(ConvexBody& body, VT& position, NT step) {
            for (int attempt = 0; attempt < 2; ++attempt) {
                VT trial = position + step * velocity_buffer +
                          NT(0.5) * step * step * acceleration_buffer;
                if (!body.isExterior(trial)) {
                    position = trial;
                    velocity_buffer += step * acceleration_buffer;
                    body.update_C(step);
                    return true;
                }
                step *= NT(0.5);
            }
            return false;
        }

        /// Metropolis acceptance test
        bool metropolisAccept(NT H0, const VT& x_prop) {
            RandomNumberGenerator &rng = settings.randomNumberGenerator;
            const NT H1 = computeHamiltonian(x_prop, velocity_buffer);
            const NT dH = H1 - H0;

            if (!std::isfinite(static_cast<double>(dH))) return false;
            if (dH <= NT(0)) return true;

            const double u = std::max(1e-300, rng.sample_urdist());
            return (-std::log(u)) > static_cast<double>(dH);
        }

        /// Adaptive step size for CONVEX optimization
        /// For convex problems: minimize reflections, maximize progress
        void adaptStepSize(NT intended, NT actual, unsigned int reflections) {
            if (!adapt_enabled) return;

            const unsigned int d = position_buffer.size();

            const NT completion = (intended > kEps()) ? (actual / intended) : NT(0);
            NT norm_dist = std::max(NT(1), actual / settings.diameter);
            const NT reflection_rate = NT(reflections) / norm_dist;

            // Dimension-adaptive thresholds
            NT good_completion = (d >= 100) ? NT(0.65) : NT(0.80);
            NT poor_completion = (d >= 100) ? NT(0.30) : NT(0.45);
            NT high_reflections = (d >= 100) ? NT(0.8) : NT(0.5);

            NT adjustment = NT(0);

            // Primary signal: completion ratio
            if (completion > good_completion && reflection_rate < high_reflections) {
                // Good progress, few reflections: increase step
                adjustment = (total_steps < 150) ? NT(0.12) : NT(0.04);
            }
            else if (completion < poor_completion || reflection_rate > high_reflections) {
                // Poor progress or too many reflections: decrease step
                adjustment = (total_steps < 150) ? -NT(0.18) : -NT(0.06);
            }
            else {
                // Middle range: gentle adjustment toward target
                const NT target = (good_completion + poor_completion) / NT(2);
                adjustment = NT(0.02) * (completion - target);
            }

            // Apply update
            const NT mult = std::exp(adjustment);
            current_step_scale = clamp(current_step_scale * mult, kMinScale(), kMaxScale());
        }
    };
};

#endif // VOLESTI_BOLTZMANN_HMC_WALK_HPP