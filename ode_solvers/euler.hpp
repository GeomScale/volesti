// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_EULER_HPP
#define ODE_SOLVERS_EULER_HPP

template <
		typename Point,
		typename NT,
		typename Polytope,
		typename func
>
struct EulerODESolver {

  typedef std::vector<Point> pts;
  typedef std::vector<Polytope*> bounds;
  typedef typename Polytope::VT VT;

  unsigned int dim;

  NT eta;
  NT t;
  VT Ar, Av;

  func F;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  // Previous state boundary point
  Point x_prev_bound;

  // Previous state boundary facet
  int prev_facet = -1;

  EulerODESolver(NT initial_time, NT step, pts initial_state, func oracle,
    bounds boundaries) :
    eta(step), t(initial_time), F(oracle), Ks(boundaries), xs(initial_state) {
      dim = xs[0].dimension();
    };


  void step(int k, bool accepted) {
    xs_prev = xs;
		t += eta;
    for (unsigned int i = 0; i < xs.size(); i++) {
      Point y = F(i, xs_prev, t);
      y = eta * y;

      if (Ks[i] == NULL) {
        xs[i] = xs_prev[i] + y;
        if (prev_facet != -1 && i > 0) {
					Ks[i-1]->compute_reflection(xs[i], x_prev_bound, prev_facet);
				}

      }
      else {

        // Find intersection (assuming a line trajectory) between x and y
        do {
          // Find line intersection between xs[i] (new position) and y
          std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs_prev[i], y, Ar, Av);

          if (pbpair.first >= 0 && pbpair.first <= 1) {
            // Advance to point on the boundary
            xs_prev[i] += (pbpair.first * 0.95) * y;

            // Update facet for reflection of derivative
            prev_facet = pbpair.second;
            x_prev_bound = xs_prev[i];

            // Reflect ray y on the boundary point y now is the reflected ray
            Ks[i]->compute_reflection(y, xs_prev[i], pbpair.second);
            // Add it to the existing (boundary) point and repeat
            xs[i] = xs_prev[i] + y;

          }
          else {
            prev_facet = -1;
            xs[i] = xs_prev[i] + y;
          }
        } while (!Ks[i]->is_in(xs[i]));

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

  void disable_adaptive() {
      // TODO Implement
  }

  void enable_adaptive() {
      // TODO Implement
  }

};

#endif
