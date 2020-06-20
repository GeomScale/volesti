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

Resource: https://en.wikipedia.org/wiki/Collocation_method

*/

#ifndef COLLOCATION_H
#define COLLOCATION_H

#include <csignal>


template <typename Point, typename NT, class Polytope, class bfunc, class func=std::function <Point(std::vector<Point>, NT)>>
class CollocationODESolver {
public:

  // Vectors of points
  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;

  // typedef from existing templates
  typedef typename Polytope::MT MT;
  typedef typename Polytope::VT VT;
  typedef std::vector<MT> MTs;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;

  unsigned int dim;

  NT eta;
  NT t, t_prev, dt;
  const NT tol = 1e-6;

  // If set to true the solver assumes linearity of the field
  // Otherwise it approximates the constant vector with Euler method
  const bool exact = false;

  // Boundary oracle method
  std::string boundary_oracle_method;

  // If set to true it enables precomputation (does not recompute A and b at every step)
  const bool precompute = true;

  bool precompute_flag = false;

  // Function oracles x'(t) = F(x, t)
  funcs Fs;

  // Basis functions
  bfunc phi, grad_phi;

  bounds Ks;

  // Contains the sub-states
  pts xs, xs_prev, zs;
  Point y;

  NT temp_grad;
  NT temp_func;

  // Basis coefficients
  ptsv as;

  // Temporal coefficients
  coeffs cs;

  // Matrices for collocation methods
  MTs As, Bs;

  // Keeps the solution to Ax = b temporarily
  MTs temps;

  VT Ar, Av;

  CollocationODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries,  coeffs c_coeffs, bfunc basis, bfunc grad_basis, std::string bmethod) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries),
     cs(c_coeffs), phi(basis), grad_phi(grad_basis), boundary_oracle_method(bmethod) {
      dim = xs[0].dimension();
      initialize_matrices();
    };

  CollocationODESolver(NT initial_time, NT step, int num_states, unsigned int dimension,
    funcs oracles, bounds boundaries,  coeffs c_coeffs,
    bfunc basis, bfunc grad_basis, std::string bmethod) :
    t(initial_time), Fs(oracles), eta(step), Ks(boundaries),  cs(c_coeffs),
    phi(basis), grad_phi(grad_basis), boundary_oracle_method(bmethod) {
      xs = pts(num_states, Point(dimension));
      initialize_matrices();
    };

  CollocationODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
     coeffs c_coeffs, bfunc basis, bfunc grad_basis, std::string bmethod) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step),
    cs(c_coeffs), phi(basis), grad_phi(grad_basis), boundary_oracle_method(bmethod) {
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
    temps = MTs(xs.size());
    as = ptsv(xs.size(), pts(order(), Point(xs[0].dimension())));
    for (unsigned int i = 0; i < xs.size(); i++) {
      // Gradient matrix is of size (order - 1) x (order - 1)
      As[i].resize(order()-1, order()-1);
      // Constants matrix is of size (order - 1) x dim
      Bs[i].resize(order()-1, xs[0].dimension());
      temps[i].resize(order()-1, xs[0].dimension());
    }

    zs = pts{Point(1)};
  }

  void step() {
    t_prev = t;
    xs_prev = xs;

    for (unsigned int ord = 0; ord < order(); ord++) {
      // Calculate t_ord
      t = t_prev + cs[ord] * eta;

      for (int i = 0; i < xs.size(); i++) {
        // a0 = F(x0, t0)
        y = Fs[i](xs, t);

        if (ord == 0) {
          as[i][0] = xs_prev[i];

          if (exact) {
            for (unsigned int j = 0; j < order() - 1; j++) {
              Bs[i].row(j) = y.getCoefficients();
            }
          }


        }
        // Construct matrix A
        else {
          if (exact) {

            if (!precompute || (precompute && !precompute_flag))  {
              for (unsigned int j = 0; j < order() - 1; j++) {
                temp_grad = grad_phi(t, t_prev, order() - j - 1, order());
                zs[0].set_coord(0, phi(t, t_prev, order() - j - 1, order()));
                temp_func = Fs[i](zs, t)[0];
                As[i](ord-1, j) = temp_grad - temp_func;
              }
            }


          }
          else {

            // Compute new derivative (inter-point)
            dt = (cs[ord] - cs[ord-1]) * eta;
            y = dt * y;


            // Do not take into account reflections
            xs[i] += y;

            // Keep grads for matrix B
            for (unsigned int j = 0; j < xs[i].dimension(); j++) {
              Bs[i].row(ord-1) = y.getCoefficients().transpose();
            }


            // Construct matrix A that contains the gradients of the basis functions
            if (!precompute || (precompute && !precompute_flag)) {
              for (unsigned int j = 0; j < order() - 1; j++) {
                As[i](ord-1, j) = grad_phi(t, t_prev, order() - j - 1, order());
              }
            }


          }
        }

        }


      }

    // Solve linear systems
    for (int i = 0; i < xs.size(); i++) {
      // temp contains solution in decreasing order of bases
      temps[i] = As[i].colPivHouseholderQr().solve(Bs[i]);

      for (int j = 0; j < order() - 1; j++) {
        // TODO Add vectorized implementation
        // as[i][order() - j - 1] += temp(j);
        for (int k = 0; k < xs[0].dimension(); k++) {
          as[i][order() - j - 1].set_coord(k, temps[i](j, k));
        }

      }
    }


    // Compute next point
    for (unsigned int i = 0; i < xs.size(); i++) {
      xs[i] = Point(xs[i].dimension());
      if (Ks[i] == NULL) {
        for (unsigned int ord = 0; ord < order(); ord++) {
          xs[i] += as[i][ord] * phi(t_prev + eta, t_prev, ord, order());
        }
      } else {
        std::tuple<NT, Point, int> result = Ks[i]->curve_intersect(t_prev, t_prev, 1.1 * eta, as[i], phi, grad_phi, boundary_oracle_method);

        std::cout << "t is " << std::get<0>(result) << std::endl;
        std::cout << "facet is " << std::get<2>(result) << std::endl;

        if (std::get<2>(result) == -1) {
          std::cout << "inside" << std::endl;
          for (unsigned int ord = 0; ord < order(); ord++) {
            xs[i] += as[i][ord] * phi(t_prev + eta, t_prev, ord, order());
          }
        }
        else {
          std::cout << "not inside" << std::endl;
          // Compute ray
          y = std::get<1>(result) - xs_prev[i];
          xs[i] = std::get<1>(result);
          // Reflect ray along facet
          Ks[i]->compute_reflection(y, xs_prev[i], std::get<2>(result));


          xs[i] += y;

          std::cout << "new point is " << xs[i].getCoefficients().transpose() << std::endl;


          while (!Ks[i]->is_in(xs[i]), 1e-6) {
            std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y, Ar, Av);

            if (pbpair.first >= 0 && pbpair.first <= 1) {
              xs[i] += (pbpair.first * 0.99) * y;
              Ks[i]->compute_reflection(y, xs[i], pbpair.second);
              xs[i] += y;
            }
            else break;

          }


        }



      }
    }

    t += eta;

    precompute_flag = true;

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

template <typename NT, class bfunc>
class RationalFunction {
  bfunc p, q;
  NT reg = 1e-6;

  RationalFunction(bfunc num, bfunc den) : p(num), q(den) {};

  NT operator()(NT t, NT t0, unsigned int j, unsigned int ord) {
    NT num = p(t, t0, j, ord);
    NT den = q(t, t0, j, ord);
    if (std::abs(den) < reg) den += reg;
    return num / den;
  }

};

template <typename NT, class bfunc>
class RationalFunctionGradient {
  bfunc p, q;
  bfunc grad_p, grad_q;
  NT reg = 1e-6;

  RationalFunctionGradient(bfunc num, bfunc grad_num, bfunc den, bfunc grad_den) :
    p(num), grad_p(grad_num), q(den), grad_q(grad_den) {};

  NT operator()(NT t, NT t0, unsigned int j, unsigned int ord) {
    NT num = p(t, t0, j, ord);
    NT grad_num = grad_p(t, t0, j, ord);
    NT den = q(t, t0, j, ord);
    NT grad_den = grad_q(t, t0, j, ord);
    if (std::abs(den * den) < reg) den += reg;
    return (grad_num  / den)  - (grad_den * num) / den;
  }

};



#endif
