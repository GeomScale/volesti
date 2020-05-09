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
*/

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

template <typename Point, typename NT>
class HamiltonianGrads {
public:
  typedef std::vector<Point> pts;
  typedef std::function<Point(pts, NT)> func;
  typedef std::vector<func> funcs;

  funcs Fs;

  HamiltonianGrads(func neg_grad_U) {
    // Oracle for position update dx / dt = v
    func grad_K = [](pts x, NT t) { return x[1]; };
    Fs.push_back(grad_K);
    // Oracle for velocity update dv / dt = - grad f
    Fs.push_back(neg_grad_U);
  }

  HamiltonianGrads() {
    // Oracle for position update dx / dt = v
    func grad_K = [](pts x, NT t) { return x[1]; };
    Fs.push_back(grad_K);
    // Oracle for velocity update dv / dt = - x
    // Uses
    func neg_grad_U = [](pts x, NT t) { return (-1.0) * x[0]; };
    Fs.push_back(neg_grad_U);
  }

  func get_grad_K() {
    return Fs[0];
  }

  func get_neg_grad_U() {
    return Fs[1];
  }

};


template <typename Point, typename NT, class Polytope>
class EulerODESolver {
public:
  typedef std::vector<Point> pts;
  typedef std::function <Point(pts, NT)> func;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;

  NT h;
  NT t;
  unsigned int counter;

  funcs Fs;
  bounds Ks;


  // Contains the sub-states
  pts xs;

  EulerODESolver(NT initial_time, NT step, pts initial_state, funcs oracles, bounds boundaries) :
    t(initial_time), xs(initial_state), Fs(oracles), counter(0), h(step), Ks(boundaries) {};

  EulerODESolver(NT initial_time, NT step, pts initial_state, funcs oracles) :
    t(initial_time), xs(initial_state), Fs(oracles), counter(0), h(step) {
      Ks = bounds(xs.size(), NULL);
    };


  void step() {
    t += h;
    NT dl = 0.95;

    for (unsigned int i = 0; i < xs.size(); i++) {
      Point y = Fs[i](xs, t);
      y = h * y;

      if (Ks[i] == NULL) {
        xs[i] = xs[i] + y;
      }
      else {
        // Find intersection (assuming a line trajectory) between x and y
        std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y);

        if (pbpair.first < 0) {
          xs[i] += (dl * pbpair.first) * y;
          Ks[i]->compute_reflection(y, xs[i], pbpair.second);
        }
        else {
          xs[i] += y;
        }
      }
      std::cout << xs[i][0] << std::endl;
    }

    counter++;
  }

  void steps(int num_steps) {
    for (int i = 0; i < num_steps; i++) step();
  }

};


#endif
