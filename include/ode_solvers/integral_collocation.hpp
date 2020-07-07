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


#ifndef COLLOCATION_HPP
#define COLLOCATION_HPP

#include "nlp_oracles/nlp_hpolyoracles.hpp"
#include "nlp_oracles/nlp_vpolyoracles.hpp"

template <
  typename Point,
  typename NT,
  class Polytope,
  class bfunc,
  class func=std::function <Point(std::vector<Point>, NT)>,
  class NontLinearOracle=MPSolveHPolyoracle<
    Polytope,
    bfunc
  >
>
class IntegralCollocationODESolver {
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

  // Function oracles x'(t) = F(x, t)
  funcs Fs;

  // Basis functions
  bfunc phi, grad_phi;

  bounds Ks;

  // Contains the sub-states
  pts xs, xs_prev;

  Point &x, &v;
  Point y;

  // Temporal coefficients
  coeffs cs;


  VT Ar, Av;

  int prev_facet = -1;
  Point prev_point;

  CollocationODESolver(NT initial_time, NT step, pts initial_state, funcs oracles,
    bounds boundaries,  coeffs c_coeffs) :
    t(initial_time), xs(initial_state), Fs(oracles), eta(step), Ks(boundaries),
     cs(c_coeffs) {
      dim = xs[0].dimension();
      initialize_matrices();
    };

  unsigned int order() const {
    return cs.size();
  }

  void initialize_matrices() {

  }

  void step() {

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
