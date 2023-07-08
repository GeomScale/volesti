// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_COLLOCATION_HPP
#define ODE_SOLVERS_COLLOCATION_HPP

#include "nlp_oracles/nlp_hpolyoracles.hpp"
#include "nlp_oracles/nlp_vpolyoracles.hpp"

template <
  typename Point,
  typename NT,
  typename Polytope,
  typename bfunc,
  typename func,
  typename NontLinearOracle=MPSolveHPolyoracle<
    Polytope,
    bfunc
  >
>
struct CollocationODESolver {


  // Vectors of points
  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;

  // typedef from existing templates
  typedef typename Polytope::MT MT;
  typedef typename Polytope::VT VT;
  typedef std::vector<MT> MTs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;

  unsigned int dim;

  NT eta;
  NT t, t_prev, dt, t_temp;
  const NT tol = 1e-3;

  // If set to true the solver assumes linearity of the field
  // Otherwise it approximates the constant vector with Euler method
  const bool exact = false;

  // If set to true it enables precomputation (does not recompute A and b at every step)
  const bool precompute = true;

  bool precompute_flag = false;

  // Function oracles x'(t) = F(x, t)
  func F;

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

  NontLinearOracle oracle;

  int prev_facet = -1;
  Point prev_point;

  CollocationODESolver(NT initial_time, NT step, pts initial_state, func oracle,
    bounds boundaries,  coeffs c_coeffs, bfunc basis, bfunc grad_basis) :
    t(initial_time), xs(initial_state), F(oracle), eta(step), Ks(boundaries),
     cs(c_coeffs), phi(basis), grad_phi(grad_basis) {
      dim = xs[0].dimension();
      initialize_matrices();
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
        y = F(i,xs, t);

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
                temp_func = F(i,zs, t)[0];
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
            Bs[i].row(ord-1) = y.getCoefficients().transpose();

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

    if (!exact) {
      for (int r = 0; r < (int) (eta / tol); r++) {
        for (unsigned int i = 0; i < xs.size(); i++) {
          xs[i] = Point(xs[i].dimension());
          for (unsigned int ord = 0; ord < order(); ord++) {
            xs[i] += as[i][ord] * phi(t_prev + eta, t_prev, ord, order());
          }
        }

        for (unsigned int i = 0; i < xs.size(); i++) {
          for (int ord = 1; ord < order(); ord++) {
            t_temp = cs[ord] * eta;
            y = F(i,xs, t_temp);
            Bs[i].row(ord-1) = y.getCoefficients().transpose();
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

      }


    }


    // Compute next point
    for (unsigned int i = 0; i < xs.size(); i++) {
      if (Ks[i] == NULL) {
        // xs[i] = Point(xs[i].dimension());
        // for (unsigned int ord = 0; ord < order(); ord++) {
        //   xs[i] += as[i][ord] * phi(t_prev + eta, t_prev, ord, order());
        // }

        if (prev_facet != -1 && i > 0) Ks[i-1]->compute_reflection(xs[i], prev_point, prev_facet);
        prev_facet = -1;


      } else {
        std::tuple<NT, Point, int> result = Ks[i]->curve_intersect(t_prev, t_prev,
            eta, as[i], phi, grad_phi, oracle);

        // Point is inside polytope
        if (std::get<2>(result) == -1 && Ks[i]->is_in(std::get<1>(result))) {
          // std::cout << "Inside" << std::endl;
          xs[i] = Point(xs[i].dimension());
          for (unsigned int ord = 0; ord < order(); ord++) {
            xs[i] += as[i][ord] * phi(t_prev + eta, t_prev, ord, order());
          }

          prev_facet = -1;

        }
        else {
          // std::cout << "outside" << std::endl;
          // Stick to the boundary
          xs[i] = 0.99 * std::get<1>(result);
          prev_point = xs[i];
          prev_facet = std::get<2>(result);


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

#endif
