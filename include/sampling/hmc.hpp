/*
GeomScale Project

Copyright (c) 2020
  Vissarion Fisikopoulos
  Apostolos Chalkis
  Elias Tsigaridas
  Marios Papachristou

Contributed and/or modified by Marios Papachristou,
as part of Google Summer of Code 2020 program.

VolEsti is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

VolEsti is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

See the file COPYING.LESSER for the text of the GNU Lesser General
Public License.  If you did not receive this file along with HeaDDaCHe,
see <http://www.gnu.org/licenses/>.

Resources:
  1. https://arxiv.org/abs/2002.04121
  2. https://arxiv.org/pdf/1801.02309
  3. http://papers.neurips.cc/paper/5801-reflection-refraction-and-hamiltonian-monte-carlo.pdf 
*/

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/erf.hpp>


#ifndef HMC_H
#define HMC_H

template <typename Point, typename NT, typename RNGType, class Polytope, class Solver>
class HMCSampler {
public:
  typedef std::vector<Point> pts;
  typedef std::function <Point(pts, NT)> func;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;

  NT eta;
  NT L, m, kappa;
  NT delta, epsilon;

  // Numerical ODE solver
  Solver *solver;

  unsigned int dim;

  // xs[0] contains position xs[1] contains velocity
  pts xs;

  // References to xs
  Point x, v;
  // Proposal points
  Point x_tilde, v_tilde;

  // Minimizer of f(x)
  Point x_min;

  // Contains K x R^d
  bounds Ks;

  // Function oracles Fs[0] contains grad_K = x
  // Fs[1] contains - grad f(x)
  funcs Fs;
  std::function<NT(Point)> f;

  boost::random::uniform_real_distribution<> urdist;
  RNGType rng;

  HMCSampler(func neg_grad_f, std::function<NT(Point)> density_exponent, Point initial, NT smoothness, NT strong_conv, NT accuracy, NT tolerance, Polytope *boundary=NULL) {
    // ODE related-stuff
    L = smoothness;
    m = strong_conv;
    kappa = L / m;
    epsilon = tolerance;
    delta = accuracy;
    eta = 1.0 / sqrt(20 * L * initial.dimension() * log(kappa / epsilon));


    // Define Kinetic and Potential Energy gradient updates
    // Kinetic energy gradient grad_K = v
    func temp_grad_K = [](pts xs, NT t) { return xs[1]; };
    Fs.push_back(temp_grad_K);
    Fs.push_back(neg_grad_f);

    // Define exp(-f(x)) where f(x) is convex
    f = density_exponent;

    // Create boundaries for K and U
    // Boundary for K is given in the constructor
    Ks.push_back(boundary);

    // Support of kinetic energy is R^d
    Ks.push_back(NULL);

    // Starting point is provided from outside
    x = initial;

    dim = initial.dimension();

    // Uniform number generator
    urdist = boost::random::uniform_real_distribution<>(0, 1);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    rng = RNGType(seed);

    solver = new Solver(0, eta, 2, initial.dimension(), Fs, Ks);


  };

  Point sample(int num_solver_steps=1) {
    // Pick a random velocity
    v = get_direction<RNGType, Point, NT>(dim, false);
    NT H = hamiltonian(x, v);
    solver->set_state(0, x);
    solver->set_state(1, v);

    // Get proposals
    solver->steps(num_solver_steps);
    x_tilde = solver->get_state(0);
    v_tilde = solver->get_state(1);


    NT H_tilde = hamiltonian(x_tilde, v_tilde);

    // Log-sum-exp trick
    NT log_prob = H - H_tilde < 0 ? H - H_tilde : 0;

    NT u_logprob = log(urdist(rng));
    if (u_logprob < log_prob) {
      x = x_tilde;
      return x;
    }
    else return x;
  }

  NT hamiltonian(Point &pos, Point &vel) {
    return f(pos) + 0.5 * vel.dot(vel);
  }


  pts samples(int num_samples, int num_solver_steps=1) {
    pts result;
    for (int i = 0; i < num_samples; i++) {
      result.push_back(sample(num_solver_steps));
    }
    return result;
  }

  void mix() {
    int steps = (int) (kappa * x.dimension() * log(1 / epsilon) * log(x.dimension() * log(kappa / epsilon)));
    samples(steps);
  }


};

#endif
