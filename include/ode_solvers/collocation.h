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

#ifndef COLLOCATION_H
#define COLLOCATION_H

template <typename Point, typename NT, class Polytope>
class CollocationSolver {
public:
  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;

  // Oracles
  typedef std::function <Point(pts, NT)> func;

  // Basis functions
  //  (point, starting point, order) -> value
  typedef std::function <NT(NT, NT, unsigned int)> bfunc;

  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;

  unsigned int dim;

  NT eta;
  NT t, t_prev;

  funcs Fs;

  // Basis functions
  bfunc phi, grad_phi;

  bounds Ks;

  // Contains the sub-states
  pts xs;
  Point y;

  // Basis coefficients
  coeffs as;

  // Temporal coefficients
  coeffs cs;


  CollocationSolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries, coeffs a_coeffs, coeffs c_coeffs, bfunc basis, bfunc grad_basis) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries),
    as(a_coeffs), cs(c_coeffs), phi(basis), grad_phi(grad_basis) {
      dim = xs[0].dimension();
    };

  CollocationSolver(NT initial_time, NT step, int num_states, unsigned int dimension,
    funcs oracles, bounds boundaries, coeffs a_coeffs, coeffs c_coeffs,
    bfunc basis, bfunc grad_basis) :
    t(initial_time), Fs(oracles), eta(step), Ks(boundaries), as(a_coeffs), cs(c_coeffs),
    phi(basis), grad_phi(grad_basis) {
      xs = pts(num_states, Point(dimension));
    };

  CollocationSolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    coeffs a_coeffs, coeffs c_coeffs, bfunc basis, bfunc grad_basis) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), as(a_coeffs),
    cs(c_coeffs), phi(basis), grad_phi(grad_basis) {
      Ks = bounds(xs.size(), NULL);
      dim = xs[0].dimension();
    };

  unsigned int order() {
    return bs.size();
  }

  void step() {


  }


  void print_state() {
    for (int j = 0; j < xs.size(); j++) {
      for (unsigned int i = 0; i < xs[j].dimension(); i++) {
        std::cout << xs[j][i] << " ";
      }
    }
    std::cout << std::endl;
  }

  void steps(int num_steps) {
    for (int i = 0; i < num_steps; i++) step();
  }

  Point get_state(int index) {
    return xs[index];
  }

  void set_state(int index, Point p) {
    xs[index] = p;
  }
};

template<typename NT>
NT poly_basis(NT t, NT t0, unsigned int j, unsigned int ord) {
  pow(t - t0, (NT) j);
}

template<typename NT>
NT poly_basis_grad(NT t, NT t0, unsigned int j, unsigned int ord) {
  ((NT) j) * pow(t - t0, (NT) (j - 1));
}

template<typename NT>
NT hermite_basis(NT t, NT t0, unsigned int ord) {

}


#endif
