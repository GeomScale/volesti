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

#ifndef BULIRSCH_STOER_H
#define BULIRSCH_STOER_H

template <typename Point, typename NT, class Polytope>
class BSODESolver {
public:
  typedef std::vector<Point> pts;
  typedef std::function <Point(pts, NT)> func;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;
  typedef std::vector<coeffs> scoeffs;
  typedef std::vector<pts> ptsv;
  typedef std::vector<ptsv> ptsm;



  unsigned int dim;
  const unsigned int MAX_TRIES = 5;

  NT eta, eta_temp;
  NT t, t_prev;
  NT tol = 1e-10;
  NT den;
  Point num, y;

  RKODESolver<Point, NT, Polytope> *solver;

  funcs Fs;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  ptsm A;
  bool flag;

  BSODESolver(NT initial_time, NT step, pts initial_state, funcs oracles, bounds boundaries) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries) {
      dim = xs[0].dimension();
      A = ptsm(MAX_TRIES, ptsv(MAX_TRIES, pts(xs.size())));
      initialize_solver();
    };


    BSODESolver(NT initial_time, NT step, int num_states, unsigned int dimension, funcs oracles, bounds boundaries) :
      t(initial_time), Fs(oracles), eta(step), Ks(boundaries) {
        for (int i = 0; i < num_states; i++) {
          xs.push_back(Point(dimension));
        }
        A = ptsm(MAX_TRIES, ptsv(MAX_TRIES, pts(num_states)));
        initialize_solver();
      };


  BSODESolver(NT initial_time, NT step, pts initial_state, funcs oracles) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step) {
      Ks = bounds(xs.size(), NULL);
      dim = xs[0].dimension();
      A = ptsm(MAX_TRIES, ptsv(MAX_TRIES, pts(xs.size())));
      initialize_solver();
    };


  void initialize_solver() {
    solver = new RKODESolver<Point, NT, Polytope>(t, eta, xs, Fs, Ks);
  }

  void step() {
    xs_prev = xs;
    eta_temp = eta;
    flag = true;

    // Use RK4 solver
    solver->xs = xs_prev;
    solver->t = t;
    solver->eta = eta_temp;
    solver->steps(1);
    A[0][0] = solver->xs;


    for (unsigned int j = 0; j < MAX_TRIES-1; j++) {
      // Reduce step size by two
      eta_temp /= 2;

      // Find solution with half stepsize and twice the num of steps
      solver->xs = xs_prev;
      solver->t = t;
      solver->eta = eta_temp;
      solver->steps(2*(1 + j));
      A[j+1][0] = solver->xs;

      // Perform Richardson extrapolation
      for (unsigned int k = 0; k <= j; k++) {
        den = 1.0 * ((4 << (k + 1)) - 1);
        for (unsigned int i = 0; i < xs.size(); i++) {
          num =  A[j+1][k][i];
          num = (1.0 * (4 << (k + 1))) * num;
          num = num - A[j][k][i];
          A[j+1][k+1][i] = (1.0 / den) * num;
        }
      }

      for (unsigned int i = 0; i < xs.size(); i++) {
        y = A[j+1][j+1][i];
        y = y - A[j][i][i];
        if (sqrt(y.dot(y)) > tol) flag = false;
      }

      if (flag) {
        for (unsigned int i = 0; i < xs.size(); i++) {
          y = A[j+1][j+1][i] - xs[i];

          if (Ks[i] == NULL) {
            xs[i] = xs[i] + y;
          }
          else {
            // Find intersection (assuming a line trajectory) between x and y
            do {
              std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y);

              if (pbpair.first < 0) {
                xs[i] += (pbpair.first * 0.99) * y;
                Ks[i]->compute_reflection(y, xs[i], pbpair.second);
              }
              else {
                xs[i] += y;
              }
            } while (!Ks[i]->is_in(xs[i]));

          }
          
        }
        break;
      }
    }

    if (!flag) {
      for (unsigned int i = 0; i < xs.size(); i++) {
        y = A[MAX_TRIES-1][MAX_TRIES-1][i] - xs[i];

        if (Ks[i] == NULL) {
          xs[i] = xs[i] + y;
        }
        else {
          // Find intersection (assuming a line trajectory) between x and y
          do {
            std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y);

            if (pbpair.first < 0) {
              xs[i] += (pbpair.first * 0.99) * y;
              Ks[i]->compute_reflection(y, xs[i], pbpair.second);
            }
            else {
              xs[i] += y;
            }
          } while (!Ks[i]->is_in(xs[i]));

        }

      }
    }

    t += eta;

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
