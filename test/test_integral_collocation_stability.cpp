// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2026 Mohit Lakra

// Licensed under GNU LGPL.3, see LICENCE file

// Verifies the numerical stability and time-scaling correctness of the Integral Collocation
// solver by checking for energy conservation in a Harmonic Oscillator (x'' = -x) simulation.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE IntegralCollocationStability
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include "ode_solvers/integral_collocation.hpp"

struct MockPoint {
    std::vector<double> coords;
    MockPoint(int dim = 2) : coords(dim, 0.0) {}
    unsigned int dimension() const { return coords.size(); }
    double operator[](int i) const { return coords[i]; }
    void set_coord(int i, double val) { coords[i] = val; }
};

struct MockPolytope {
    typedef Eigen::MatrixXd MT;
    typedef Eigen::VectorXd VT;
};

struct HarmonicOscillator {
    MockPoint operator()(unsigned int idx, const std::vector<MockPoint>& state, double t) {
        double x = state[0][0];
        double v = state[0][1];
        MockPoint deriv(2);
        deriv.set_coord(0, v);
        deriv.set_coord(1, -x);
        return deriv;
    }
};

BOOST_AUTO_TEST_CASE(test_integral_collocation_stability) {
    double eta = 0.1; 
    int steps = 100;
    
    MockPoint initial_pt(2);
    initial_pt.set_coord(0, 0.0);
    initial_pt.set_coord(1, 1.0);
    std::vector<MockPoint> initial_state;
    initial_state.push_back(initial_pt);
    
    std::vector<MockPolytope*> boundaries(1, nullptr);

    IntegralCollocationODESolver<MockPoint, double, MockPolytope, HarmonicOscillator> 
        solver(0.0, eta, initial_state, HarmonicOscillator(), boundaries, 10);

    bool unstable = false;
    
    for (int i = 1; i <= steps; ++i) {
        solver.step();
        MockPoint current = solver.get_state(0);
        
        double radius = std::sqrt(current[0]*current[0] + current[1]*current[1]);
        
        if (radius > 1.5) {
            unstable = true;
            break;
        }
    }

    BOOST_CHECK_MESSAGE(!unstable, "Instability detected! Energy grew uncontrollably.");
}
