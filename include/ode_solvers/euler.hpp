// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef EULER_HPP
#define EULER_HPP

template <typename Point, typename NT, class Polytope, class func=std::function <Point(std::vector<Point>, NT)>>
class EulerODESolver {
public:
  typedef std::vector<Point> pts;
  typedef std::vector<func> funcs;
  typedef std::vector<Polytope*> bounds;
  typedef typename Polytope::VT VT;

  unsigned int dim;

  NT eta;
  NT t;
  VT Ar, Av;

  funcs Fs;
  bounds Ks;

  // Contains the sub-states
  pts xs;
  pts xs_prev;

  EulerODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries) {
      dim = xs[0].dimension();
    };


    EulerODESolver(NT initial_time, NT step, int num_states,
      unsigned int dimension, funcs oracles, bounds boundaries) :
      t(initial_time), Fs(oracles), eta(step), Ks(boundaries) {
        for (int i = 0; i < num_states; i++) {
          xs.push_back(Point(dimension));
        }
      };


  EulerODESolver(NT initial_time, NT step, pts initial_state, funcs oracles) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step) {
      Ks = bounds(xs.size(), NULL);
      dim = xs[0].dimension();
    };


  void step() {
    xs_prev = xs;
    t += eta;

    bool flag;

    for (unsigned int i = 0; i < xs.size(); i++) {
      Point y = Fs[i](xs_prev, t);
      y = eta * y;

      if (Ks[i] == NULL) {
        xs[i] = xs[i] + y;
      }
      else {
        flag = false;

        // Find intersection (assuming a line trajectory) between x and y
        do {
          // Find line intersection between xs[i] (new position) and y
          std::pair<NT, int> pbpair = Ks[i]->line_positive_intersect(xs[i], y, Ar, Av);

          // If point is outside it would yield a negative param
          if (pbpair.first >= 0 && pbpair.first <= 1) {
            // Advance to point on the boundary
            xs[i] += (pbpair.first * 0.99) * y;
            // Reflect ray y on the boundary point y now is the reflected ray
            Ks[i]->compute_reflection(y, xs[i], pbpair.second);
            // Add it to the existing (boundary) point and repeat
            xs[i] += y;
          }
          else {
            if (flag) break;
            xs[i] += y;
            flag = true;
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
