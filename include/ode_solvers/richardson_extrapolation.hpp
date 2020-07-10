// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RICHARDSON_EXTRAPOLATION_HPP
#define RICHARDSON_EXTRAPOLATION_HPP


template <typename Point, typename NT, class Polytope, class func=std::function <Point(std::vector<Point>&, NT&)>>
class RichardsonExtrapolationODESolver {
public:
  typedef std::vector<Point> pts;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef std::vector<NT> coeffs;
  typedef std::vector<coeffs> scoeffs;
  typedef std::vector<pts> ptsv;
  typedef std::vector<ptsv> ptsm;

  typedef typename Polytope::VT VT;

  unsigned int dim;
  const unsigned int MAX_TRIES = 5;

  NT eta, eta_temp;
  NT t, t_prev;
  NT tol = 1e-7;
  NT error = NT(-1);
  NT den;
  Point num, y;
  VT Ar, Av;

  RKODESolver<Point, NT, Polytope> *solver;

  funcs Fs;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  ptsm A;
  bool flag;

  // Previous state boundary point
  Point x_prev_bound;

  // Previous state boundary facet
  int prev_facet = -1;

  RichardsonExtrapolationODESolver(NT initial_time, NT step, pts initial_state,
    funcs oracles, bounds boundaries) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries) {
      dim = xs[0].dimension();
      A = ptsm(MAX_TRIES+1, ptsv(MAX_TRIES+1, pts(xs.size())));
      initialize_solver();
    };

  void initialize_solver() {
    solver = new RKODESolver<Point, NT, Polytope>(t, eta, xs, Fs, bounds{NULL});
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
    A[1][1] = solver->xs;


    for (unsigned int j = 1; j <= MAX_TRIES-1; j++) {
      // Reduce step size by two
      eta_temp /= 2;

      // Find solution with half stepsize and twice the num of steps
      solver->xs = xs_prev;
      solver->t = t;
      solver->eta = eta_temp;
      solver->steps(2*j);
      A[j+1][1] = solver->xs;

      // Perform Richardson extrapolation
      for (unsigned int k = 1; k <= j; k++) {
        den = 1.0 * ((4 << k) - 1);
        for (unsigned int i = 0; i < xs.size(); i++) {
          num = (1.0 * (4 << k)) * A[j+1][k][i];
          num = num - A[j][k][i];
          A[j+1][k+1][i] = (1 / den) * num;
        }
      }

      error = NT(-1);

      for (unsigned int i = 0; i < xs.size(); i++) {
        y = A[j+1][j+1][i] - A[j][j][i];
        if (sqrt(y.dot(y)) > error) error = sqrt(y.dot(y));
      }

      if (error < tol || j == MAX_TRIES - 1) {

        for (unsigned int i = 0; i < xs.size(); i++) {
          y = A[j+1][j+1][i] - xs[i];

          if (Ks[i] == NULL) {
            xs[i] = xs[i] + y;
            if (prev_facet != -1) Ks[i]->compute_reflection(xs[i], x_prev_bound, prev_facet);
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
        break;
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
