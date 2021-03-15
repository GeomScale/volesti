// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_RUNGE_KUTTA_H
#define ODE_SOLVERS_RUNGE_KUTTA_H


template <
		typename Point,
		typename NT,
		typename Polytope,
		typename func
>
struct RKODESolver {

  typedef std::vector<Point> pts;
  typedef std::vector<pts> ptsv;

  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;
  typedef std::vector<coeffs> scoeffs;

  typedef typename Polytope::VT VT;

  unsigned int dim;

  NT eta;
  NT t, t_prev;

  VT Ar, Av;

  func F;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  ptsv ks;
  Point y;

  // Previous state boundary point
  Point x_prev_bound;

  // Previous state boundary facet
  int prev_facet = -1;

  scoeffs as;
  coeffs cs, bs;


  RKODESolver(NT initial_time, NT step, pts initial_state, func oracle,
    bounds boundaries) :
    eta(step), t(initial_time), F(oracle), Ks(boundaries), xs(initial_state) {
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


  // RKODESolver(NT initial_time, NT step, pts initial_state, func oracle,
  //   bounds boundaries, scoeffs a_coeffs, coeffs b_coeffs, coeffs c_coeffs) :
  //   t(initial_time), xs(initial_state), F(oracle), eta(step), Ks(boundaries),
  //   as(a_coeffs), bs(b_coeffs), cs(c_coeffs) {
  //     dim = xs[0].dimension();
  //   };

  unsigned int order() {
    return bs.size();
  }

  void step(int k, bool accepted) {
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
        y = F(i,ks[ord], t);
        ks[ord][i] = y;
        y = (eta * bs[i]) * y;

        if (Ks[i] == NULL) {
          xs[i] = xs[i] + y;
          if (prev_facet != -1 && i > 0)
						Ks[i-1]->compute_reflection(xs[i], x_prev_bound, prev_facet);
          prev_facet = -1;
        }
        else {

          // Find intersection (assuming a line trajectory) between x and y
          do {
            // Find line intersection between xs[i] (new position) and y
            std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y, Ar, Av);
            // If point is outside it would yield a negative param
            if (pbpair.first >= 0 && pbpair.first <= 1) {
              // Advance to point on the boundary
              xs[i] += (pbpair.first * 0.99) * y;

              // Update facet for reflection of derivative
              prev_facet = pbpair.second;
              x_prev_bound = xs[i];

              // Reflect ray y on the boundary point y now is the reflected ray
              Ks[i]->compute_reflection(y, xs[i], pbpair.second);
              // Add it to the existing (boundary) point and repeat
              xs[i] += y;

            }
            else {
              prev_facet = -1;
              xs[i] += y;
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

  void steps(int num_steps, bool accepted) {
    for (int i = 0; i < num_steps; i++) step(i, accepted);
  }

  Point get_state(int index) {
    return xs[index];
  }

  void set_state(int index, Point p) {
    xs[index] = p;
  }

};


#endif
