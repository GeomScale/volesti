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

Resource: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

*/

#ifndef RK_H
#define RK_H


template <typename Point, typename NT, class Polytope, class func=std::function <Point(std::vector<Point>, NT)>>
class RKODESolver {
public:
  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;
  typedef std::vector<coeffs> scoeffs;

  typedef typename Polytope::VT VT;

  unsigned int dim;

  NT eta;
  NT t, t_prev;

  VT Ar, Av;

  funcs Fs;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  ptsv ks;
  Point y;

  scoeffs as;
  coeffs cs, bs;


  RKODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries) {
      dim = xs[0].dimension();
      // If no coefficients are given the RK4 method is assumed
      cs = coeffs{0, 0.5, 0.5, 1};
      bs = coeffs{1.0/6, 1.0/3, 1.0/3, 1.0/6};
      as = scoeffs{
        coeffs{},
        coeffs{0.5},
        coeffs{0, 0.5},
        coeffs{0, 0, 1.0}
      };
    };


  RKODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries, scoeffs a_coeffs, coeffs b_coeffs, coeffs c_coeffs) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries),
    as(a_coeffs), bs(b_coeffs), cs(c_coeffs) {
      dim = xs[0].dimension();
    };




  RKODESolver(NT initial_time, NT step, int num_states, unsigned int dimension,
    funcs oracles, bounds boundaries, scoeffs a_coeffs, coeffs b_coeffs, coeffs c_coeffs) :
    t(initial_time), Fs(oracles), eta(step), Ks(boundaries), as(a_coeffs),
    bs(b_coeffs), cs(c_coeffs) {
      xs = pts(num_states, Point(dimension));
    };

  RKODESolver(NT initial_time, NT step, int num_states, unsigned int dimension,
    funcs oracles, bounds boundaries) :
    t(initial_time), Fs(oracles), eta(step), Ks(boundaries) {
      xs = pts(num_states, Point(dimension));
      // If no coefficients are given the RK4 method is assumed
      cs = coeffs{0, 0.5, 0.5, 1};
      bs = coeffs{1.0/6, 1.0/3, 1.0/3, 1.0/6};
      as = scoeffs{
        coeffs{},
        coeffs{0.5},
        coeffs{0, 0.5},
        coeffs{0, 0, 1.0}
      };
    };


  RKODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    scoeffs a_coeffs, coeffs b_coeffs, coeffs c_coeffs) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), as(a_coeffs),
    bs(b_coeffs), cs(c_coeffs) {
      Ks = bounds(xs.size(), NULL);
      dim = xs[0].dimension();
    };

  RKODESolver(NT initial_time, NT step, pts initial_state, funcs oracles) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step) {
      Ks = bounds(xs.size(), NULL);
      dim = xs[0].dimension();
      // If no coefficients are given the RK4 method is assumed
      cs = coeffs{0, 0.5, 0.5, 1};
      bs = coeffs{1.0/6, 1.0/3, 1.0/3, 1.0/6};
      as = scoeffs{
        coeffs{},
        coeffs{0.5},
        coeffs{0, 0.5},
        coeffs{0, 0, 1.0}
      };
    };


  unsigned int order() {
    return bs.size();
  }

  void step() {
    ks = ptsv(order(), xs);
    t_prev = t;

    for (unsigned int ord = 0; ord < order(); ord++) {
      // Initialize t to previous
      t = t_prev + cs[ord] * eta;

      // Initialize ks to previous solution we use
      // Initialize argument
      for (int j = 0; j < ord; j++) {
        for (int r = 0; r < xs.size(); r++) {
          y =  ks[j][r];
          y = (eta * as[ord][j]) * y;
          ks[ord][r] = ks[ord][r] + y;
        }
      }

      for (unsigned int i = 0; i < xs.size(); i++) {
        // Calculate k_i s
        y = Fs[i](ks[ord], t);
        ks[ord][i] = y;
        y = (eta * bs[i]) * y;

        if (Ks[i] == NULL) {
          xs[i] = xs[i] + y;

        }
        else {
          // Find intersection (assuming a line trajectory) between x and y
          do {
            std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y, Ar, Av);

            if (pbpair.first < 0) {
              xs[i] += (pbpair.first * 0.99) * y;
              Ks[i]->compute_reflection(y, xs[i], pbpair.second);
            }
            else {
              xs[i] += y;
              break;
            }
          } while (!Ks[i]->is_in(xs[i]));
        }
      }
    }

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


#endif
