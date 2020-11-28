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

  std::vector<VT> Ar, Av;
  std::vector<NT> lambda_prev;

  NT eta;
  NT eta0;
  NT t;
  NT dl = 0.95;

  func F;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  std::pair<NT, int> pbpair;

  unsigned long long num_reflections = 0;
  unsigned long long num_steps = 0;

  bool adaptive = true;

  LeapfrogODESolver(NT initial_time, NT step, pts initial_state, func oracle, bounds boundaries, bool adaptive_=true) :
  eta(step), eta0(step), t(initial_time), F(oracle), Ks(boundaries), xs(initial_state), adaptive(adaptive_) {
    dim = xs[0].dimension();
    initialize();
  };

  void initialize() {
      for (unsigned int i = 0; i < xs.size(); i++) {
          VT ar, av;
          if (Ks[i] != NULL) {
              ar.setZero(Ks[i]->num_of_hyperplanes());
              av.setZero(Ks[i]->num_of_hyperplanes());
          }
          Ar.push_back(ar);
          Av.push_back(av);
          lambda_prev.push_back(NT(0));
      }
      step();
  }

  void step() {
    num_steps++;

    if (adaptive) eta = (eta0 * num_steps) / (num_steps + num_reflections);

    xs_prev = xs;
    unsigned int x_index, v_index, it;
    t += eta;
    for (unsigned int i = 1; i < xs.size(); i += 2) {
        pbpair.second  = -1;
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
        unsigned int it = 0;
        do {
          pbpair = Ks[x_index]->line_positive_intersect(xs_prev[x_index], y, Ar[x_index], Av[x_index]);

          if (pbpair.first >= 0 && pbpair.first <= 1) {
            num_reflections++;

            xs_prev[x_index] += (pbpair.first * dl) * y;
            lambda_prev[x_index] = dl * pbpair.first;
            Ks[x_index]->compute_reflection(y, xs_prev[x_index], pbpair.second);

            y = (1 - dl * pbpair.first) * y;
            xs[x_index] = xs_prev[x_index] + y;

            // Reflect velocity
            Ks[x_index]->compute_reflection(xs[v_index], xs[x_index], pbpair.second);
            xs[v_index] = (1 - dl * pbpair.first) * xs[v_index];
          }
          else {
            xs[x_index] = xs_prev[x_index] + y;
            break;
          }
      } while (!Ks[x_index]->is_in(xs[x_index]));
      }

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
