// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef LEAPFROG_HPP
#define LEAPFROG_HPP

template <
typename Point,
typename NT,
typename Polytope,
typename func
>
struct LeapfrogODESolver {

  typedef std::vector<Point> pts;

  typedef std::vector<Polytope*> bounds;
  typedef typename Polytope::VT VT;

  unsigned int dim;

  VT Ar, Av;

  NT eta;
  NT eta0;
  NT t;

  func F;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;
  pts xs_prev_backup;

  std::pair<NT, int> pbpair;

  unsigned long long num_reflections = 0;
  unsigned long long num_steps = 0;

  bool adaptive = true;

  LeapfrogODESolver(NT initial_time, NT step, pts initial_state, func oracle, bounds boundaries, bool adaptive_=true) :
  eta(step), eta0(step), t(initial_time), F(oracle), Ks(boundaries), xs(initial_state), adaptive(adaptive_) {
    dim = xs[0].dimension();
  };

  void step() {
    num_steps++;

    if (adaptive) eta = (eta0 * num_steps) / (num_steps + num_reflections);

    xs_prev = xs;
    xs_prev_backup = xs;

    unsigned int x_index, v_index, it;
    t += eta;

    for (unsigned int i = 1; i < xs.size(); i += 2) {
      x_index = i - 1;
      v_index = i;

      // v' <- v + eta / 2 F(x)
      Point z = F(v_index, xs_prev, t);
      z = (eta / 2) * z;
      xs[v_index] = xs[v_index] + z;

      // x <- x + eta v'
      Point y = xs[v_index];
      y = (eta) * y;

      if (Ks[x_index] == NULL) {
        xs[x_index] = xs_prev[x_index] + y;
      }
      else {
        // Find intersection (assuming a line trajectory) between x and y
        it = 0;
        do {
          if (it > 100 * dim) {
            xs = xs_prev_backup;
            return;
          }
          pbpair = Ks[x_index]->line_positive_intersect(xs_prev[x_index], y, Ar, Av);

          if (pbpair.first >= 0 && pbpair.first <= 1) {
            it++;
            xs_prev[x_index] += (pbpair.first * 0.95) * y;
            Ks[x_index]->compute_reflection(y, xs_prev[x_index], pbpair.second);
            xs[x_index] = xs_prev[x_index] + y;

            // Reflect velocity
            Ks[x_index]->compute_reflection(xs[v_index], xs[x_index], pbpair.second);
          }
          else {
            xs[x_index] = xs_prev[x_index] + y;
          }
        } while (!Ks[x_index]->is_in(xs[x_index]));
      }

      num_reflections += it;

      // tilde v <- v + eta / 2 F(tilde x)
      z = F(v_index, xs, t);
      z = (eta / 2) * z;
      xs[v_index] = xs[v_index] + z;

    }

  }

  void print_state() {
    for (int j = 0; j < xs.size(); j ++) {
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
