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

template <typename Point, typename NT, class Polytope, class bfunc>
class CollocationODESolver {
public:
  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;

  // Coefficient matrices
  typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
  typedef std::vector<MT> MTs;

  // Oracles
  typedef std::function <Point(pts, NT)> func;


  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;

  unsigned int dim;

  NT eta;
  NT t, t_prev, dt;

  funcs Fs;

  // Basis functions
  bfunc phi, grad_phi;

  bounds Ks;

  // Contains the sub-states
  pts xs, xs_prev;
  Point y;

  // Basis coefficients
  ptsv as;

  // Temporal coefficients
  coeffs cs;

  // Matrices for collocation methods
  MTs As, Bs;

  // Keeps the solution to Ax = b temporarily
  MT temp;

  CollocationODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries,  coeffs c_coeffs, bfunc basis, bfunc grad_basis) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries),
     cs(c_coeffs), phi(basis), grad_phi(grad_basis) {
      dim = xs[0].dimension();
      initialize_matrices();
    };

  CollocationODESolver(NT initial_time, NT step, int num_states, unsigned int dimension,
    funcs oracles, bounds boundaries,  coeffs c_coeffs,
    bfunc basis, bfunc grad_basis) :
    t(initial_time), Fs(oracles), eta(step), Ks(boundaries),  cs(c_coeffs),
    phi(basis), grad_phi(grad_basis) {
      xs = pts(num_states, Point(dimension));
      initialize_matrices();
    };

  CollocationODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
     coeffs c_coeffs, bfunc basis, bfunc grad_basis) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step),
    cs(c_coeffs), phi(basis), grad_phi(grad_basis) {
      Ks = bounds(xs.size(), NULL);
      initialize_matrices();
      dim = xs[0].dimension();
    };

  unsigned int order() const {
    return cs.size();
  }

  void initialize_matrices() {
    As = MTs(xs.size());
    Bs = MTs(xs.size());
    as = ptsv(xs.size(), pts(order(), Point(xs[0].dimension())));
    for (unsigned int i = 0; i < xs.size(); i++) {
      // Gradient matrix is of size (order - 1) x (order - 1)
      As[i].resize(order()-1, order()-1);
      // Constants matrix is of size (order - 1) x dim
      Bs[i].resize(order()-1, xs[0].dimension());
    }
    temp.resize(order()-1, xs[0].dimension());
  }

  void step() {
    t_prev = t;

    for (unsigned int ord = 0; ord < order(); ord++) {
      // Calculate t_ord
      t = t_prev + cs[ord] * eta;

      for (int i = 0; i < xs.size(); i++) {
        y = Fs[i](xs, t);

        // a0 = F(x0, t0)
        if (ord == 0) as[i][0] = y;
        else {
          // Construct matrix b
          dt = (cs[ord] - cs[ord-1]) * eta;

          // Compute new derivative (inter-point)
          y = dt * y;

          // Do not take into account reflections
          xs[i] += y;

          // Construct matrix A that contains the gradients of the basis functions
          for (unsigned int j = 0; j < order() - 1; j++) {
            As[i](ord-1, j) = grad_phi(t, t_prev, order() - j - 1, order());
          }

          // Keep grads for matrix B
          for (unsigned int j = 0; j < xs[i].dimension(); j++) {
            Bs[i](ord-1, j) = y[j];
          }

        }
      }
    }


    // Solve linear systems
    for (int i = 0; i < xs.size(); i++) {
      temp = As[i].colPivHouseholderQr().solve(Bs[i]);
      for (int j = 1; j < order(); j++) {
        for (int k = 0; k < xs[0].dimension(); k++) {
          as[i][j].set_coord(k, temp(j-1, k));
        }
      }
    }

    // Compute next point
    for (unsigned int i = 0; i < xs.size(); i++) {
      if (Ks[i] == NULL) {
        for (unsigned int ord = 0; ord < order(); ord++) {
          xs[i] += as[i][ord] * phi(t_prev + eta, t_prev, ord, order());
        }
      } else {
        // TODO implement
        throw true;

      }
    }

    print_state();
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
  return pow(t - t0, (NT) j);
}

template<typename NT>
NT poly_basis_grad(NT t, NT t0, unsigned int j, unsigned int ord) {
  return ((NT) j) * pow(t - t0, (NT) (j - 1));
}


#endif
