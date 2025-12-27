// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

// Refers to the integral collocation method with Lagrange Polynomials
// from Lee, Yin Tat, Zhao Song, and Santosh S. Vempala.
//"Algorithmic theory of ODEs and sampling from well-conditioned
// logconcave densities." arXiv preprint arXiv:1812.06243 (2018).


#ifndef ODE_SOLVERS_INTEGRAL_COLLOCATION_HPP
#define ODE_SOLVERS_INTEGRAL_COLLOCATION_HPP

#include <cmath>
#include <limits>

#include "nlp_oracles/nlp_hpolyoracles.hpp"
#include "nlp_oracles/nlp_vpolyoracles.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/math/special_functions/chebyshev.hpp"
#include "boost/math/special_functions/chebyshev_transform.hpp"
#include "root_finders/mp_solve_wrapper.hpp"

template <
    typename Point,
    typename NT,
    typename Polytope,
    typename func
>
struct IntegralCollocationODESolver {


  // Vectors of points
  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;

  // typedef from existing templates
  typedef typename Polytope::MT MT;
  typedef typename Polytope::VT VT;
  typedef std::vector<MT> MTs;

  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;
  typedef boost::numeric::ublas::vector<NT> boost_vector;
  typedef boost::math::chebyshev_transform<NT> chebyshev_transform_boost;

  unsigned int dim;

  NT eta;
  NT t, t_prev, dt, temp_node, a, b;
  const NT tol = 1e-6;

  // Function oracles x'(t) = F(x, t)
  func F;
  bounds Ks;

  // Contains the sub-states
  pts xs, xs_prev, X_temp;
  Point y;

  // Temporal coefficients
  coeffs cs;

  VT Ar, Av, X_op, nodes;

  MT A_phi, X0, X, X_prev, F_op;

  unsigned int _order;

  LagrangePolynomial<NT, VT> lagrange_poly;

  int prev_facet = -1;
  Point prev_point;

  Point evaluate_polynomial(const pts &coefficients, NT value) const {
    if (coefficients.empty()) {
      return Point(dim);
    }

    Point result(coefficients[0].dimension());
    NT current_power = NT(1);

    for (const Point &coeff : coefficients) {
      result += coeff * current_power;
      current_power *= value;
    }

    return result;
  }

  IntegralCollocationODESolver(NT initial_time, NT step, pts initial_state,
    func oracle, bounds boundaries, unsigned int order_=4) :
    t(initial_time), xs(initial_state), X_temp(initial_state), F(oracle), eta(step), Ks(boundaries),
    _order(order_) {
      dim = xs[0].dimension();
      initialize_matrices();
    };

  unsigned int order() const {
    return _order;
  }

  void initialize_matrices() {

    A_phi.resize(order(), order());
    nodes.resize(order());

    std::vector<NT> temp;

    for (unsigned int j = 0; j < order(); j++) {
      nodes(j) = cos((j+0.5) * M_PI / order());
    }

    lagrange_poly.set_nodes(nodes);

    // Calculate integrals of basis functions based on the Discrete Chebyshev Transform
    for (unsigned int i = 0; i < order(); i++) {

      lagrange_poly.set_basis((int) i);

      for (unsigned int j = 0; j <= i; j++) {
        if (nodes(j) < NT(0)) {
          a = nodes(j);
          b = NT(0);
        } else {
          a = NT(0);
          b = nodes(j);
        }

        chebyshev_transform_boost transform(lagrange_poly, a, b);
        A_phi(i, j) =  NT(transform.integrate());
        A_phi(j, i) = A_phi(i, j);
      }
    }

    #ifdef VOLESTI_DEBUG
      std::cout << "A_phi" << std::endl;
      std::cout << A_phi << std::endl;
    #endif

    X.resize(xs.size() * dim, order());
    X0.resize(xs.size() * dim, order());
    X_prev.resize(xs.size() * dim, order());
    X_op.resize(xs.size() * dim);

    F_op.resize(xs.size() * dim, order());

    lagrange_poly.set_basis(-1);
  }

  void initialize_fixed_point() {
    for (unsigned int ord = 0; ord < order(); ord++) {
      for (unsigned int i = 0; i < xs.size(); i++) {
        for (unsigned int j = i * dim; j < (i + 1) * dim; j++) {
          X0(j, ord) = xs_prev[i][j % dim];
        }
      }
    }
  }

  void step() {
    xs_prev = xs;
    initialize_fixed_point();

    std::vector<chebyshev_transform_boost> transforms;

    X = X0;
    X_prev = 100 * X0;
    NT err;

    do {
      for (unsigned int ord = 0; ord < order(); ord++) {
        for (unsigned int i = 0; i < xs.size(); i++) {
          for (unsigned int j = i * dim; j < (i + 1) * dim; j++) {
            X_temp[i].set_coord(j % dim, X(j, ord));
          }
        }

        for (unsigned int i = 0; i < xs.size(); i++) {
          // std::cout << "pre y" << std::endl;
          temp_node = nodes(ord) * eta;
          y = F(i, X_temp, temp_node);

          for (int j = i * dim; j < (i + 1) * dim; j++) {
            F_op(j, ord) = y[j % dim];
          }

        }
      }

      X = X0 + F_op * A_phi;

      X_prev = X;

      err = sqrt((X - X_prev).squaredNorm());

    } while (err > 1e-10);

    X_op = X0.col(0);

    unsigned int max_transform_coeffs_length = 0;


    for (unsigned int i = 0; i < xs.size(); i++) {
      for (unsigned int j = i * dim; j < (i + 1) * dim; j++) {
        lagrange_poly.set_coeffs(F_op.row(j).transpose());
        chebyshev_transform_boost transform(lagrange_poly, 0, eta, 1e-5, 5);
        transforms.push_back(transform);
        // Keep max transform length for zero-padding
        if (max_transform_coeffs_length < transform.coefficients().size()) {
          max_transform_coeffs_length = transform.coefficients().size();
        }

        X_op(j) += NT(transform.integrate());
      }
    }


    for (unsigned int i = 0; i < xs.size(); i++) {
      if (Ks[i] == NULL) {
        for (unsigned int j = i * dim; j < (i + 1) * dim; j++) {
          xs[i].set_coord(j % dim, X_op(j));
        }
      }
      else {
        if (max_transform_coeffs_length == 0) {
          for (unsigned int j = i * dim; j < (i + 1) * dim; j++) {
            xs[i].set_coord(j % dim, X_op(j));
          }
          continue;
        }

        // Build a polynomial representation of the trajectory on [0, eta].
        pts chebyshev_coeffs(max_transform_coeffs_length, Point(dim));
        for (unsigned int coord = 0; coord < dim; coord++) {
          unsigned int transform_index = i * dim + coord;
          std::vector<NT> temp_coeffs = transforms[transform_index].coefficients();
          for (unsigned int k = 0; k < temp_coeffs.size(); k++) {
            chebyshev_coeffs[k].set_coord(coord, temp_coeffs[k]);
          }
        }

        unsigned int doubled_length = 2 * (max_transform_coeffs_length - 1) + 1;
        pts trajectory_coeffs(doubled_length, Point(dim));
        degree_doubling_chebyshev(chebyshev_coeffs, trajectory_coeffs);

        trajectory_coeffs.insert(trajectory_coeffs.begin(), Point(dim));
        for (unsigned int k = 1; k < trajectory_coeffs.size(); k++) {
          trajectory_coeffs[k] = (NT(1) / NT(k)) * trajectory_coeffs[k];
        }
        trajectory_coeffs[0] = xs_prev[i];

        Point candidate_point = evaluate_polynomial(trajectory_coeffs, eta);

        if (Ks[i]->is_in(candidate_point, tol) == -1) {
          xs[i] = candidate_point;
          prev_facet = -1;
          continue;
        }

        MT A_bound = Ks[i]->get_mat();
        VT b_bound = Ks[i]->get_vec();
        NT best_time = eta;
        Point best_point = candidate_point;
        int best_facet = -1;
        const NT boundary_tol = NT(1e-5);
        const NT root_tol = NT(1e-9);

        // Project the polynomial on each facet normal and locate the earliest hit.
        for (int facet = 0; facet < A_bound.rows(); facet++) {
          std::vector<NT> projected(trajectory_coeffs.size(), NT(0));
          for (unsigned int deg = 0; deg < trajectory_coeffs.size(); deg++) {
            projected[deg] = A_bound.row(facet) * trajectory_coeffs[deg].getCoefficients();
          }
          projected[0] -= b_bound(facet);

          auto roots = mpsolve<NT>(projected);
          for (auto &root : roots) {
            NT candidate_time = root.first;
            if (candidate_time < root_tol || candidate_time > eta + tol) continue;

            Point intersection = evaluate_polynomial(trajectory_coeffs, candidate_time);
            NT facet_value = A_bound.row(facet) * intersection.getCoefficients() - b_bound(facet);

            if (std::abs(facet_value) > boundary_tol) continue;
            if (Ks[i]->is_in(intersection, boundary_tol) == 0) continue;

            if (candidate_time < best_time) {
              best_time = candidate_time;
              best_point = intersection;
              best_facet = facet;
            }
          }
        }

        if (best_facet != -1) {
          Point offset = best_point - xs_prev[i];
          xs[i] = xs_prev[i] + NT(0.99) * offset;
          prev_point = best_point;
          prev_facet = best_facet;
        } else {
          Point direction = candidate_point - xs_prev[i];
          NT direction_norm = direction.getCoefficients().squaredNorm();
          if (direction_norm < tol) {
            xs[i] = xs_prev[i];
          } else {
            std::pair<NT, NT> line_hit = Ks[i]->line_intersect(xs_prev[i], direction);
            NT lambda = line_hit.first;
            if (lambda < std::numeric_limits<NT>::max()) {
              xs[i] = xs_prev[i] + (direction * (NT(0.99) * lambda));
            } else {
              xs[i] = candidate_point;
            }
          }
          prev_facet = -1;
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
