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
class EulerODESolver {
public:
  typedef std::function <Point(Point, NT)> funcs;
  // Polytope *K;
  NT h;
  NT t;
  unsigned int T;
  unsigned int counter;

  funcs F;

  // Contains the sub-states
  Point x;

  EulerODESolver(NT initial_time, unsigned int total_steps, NT step, Point initial_state, funcs oracles) :
    t(initial_time), T(total_steps), x(initial_state), F(oracles), counter(0), h(step) {};

  void step() {
    t += h;
    Point y = F(x, t);
    y = h * y;
    x = x + y;

    std::cout << x[0] << std::endl;

    // TODO add boundary reflections
    counter++;
  }

  void steps() {
    while (counter < T) step();
  }

};


#endif
