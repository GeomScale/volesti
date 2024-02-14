// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_ORACLE_FUNCTORS_HPP
#define ODE_SOLVERS_ORACLE_FUNCTORS_HPP

struct OptimizationFunctor {
    template <
        typename NT,
        typename Functor,
        typename GradFunctor
    >
    struct parameters {
        NT T; // Temperature
        unsigned int dim; // Dimension
        Functor f;
        GradFunctor neg_grad_f;
        NT L;
        NT m;
        NT kappa;
        unsigned int order;

        parameters(
            NT T_,
            unsigned int dim_,
            Functor f_,
            GradFunctor neg_grad_f_) :
            T(T_),
            dim(dim_),
            f(f_),
            neg_grad_f(neg_grad_f_),
            L(1),
            m(1),
            kappa(1),
            order(2)
        {};

        void update_temperature(NT k, NT l) {
            T = T * pow(1.0 + 1.0 / pow(dim, k), l);
        }
    };

    template
    <
        typename Point,
        typename Functor,
        typename GradFunctor
    >
    struct GradientFunctor {
        typedef typename Point::FT NT;
        typedef std::vector<Point> pts;

        parameters<NT, Functor, GradFunctor> &params;

        GradientFunctor(parameters<NT, Functor, GradFunctor> &params_) : params(params_) {};

        // The index i represents the state vector index
        Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
          if (i == params.order - 1) {
            return params.neg_grad_f(i, xs, t) * params.T; // returns - a*x
          } else {
            return xs[i + 1]; // returns derivative
          }
        }
      };

    template
    <
        typename Point,
        typename Functor,
        typename GradFunctor
    >
    struct FunctionFunctor {
        typedef typename Point::FT NT;
        parameters<NT, Functor, GradFunctor> &params;

        FunctionFunctor(parameters<NT, Functor, GradFunctor> &params_) : params(params_) {};

        NT operator() (Point const& x) const {
            return params.f(x) * params.T;
        }
    };
};

struct IsotropicQuadraticFunctor {

    // Holds function oracle and gradient oracle for the function 1/2 a ||x||^2
    template <
        typename NT
    >
    struct parameters {
        NT alpha;
        unsigned int order;
        NT L; // Lipschitz constant of gradient
        NT m; // Strong-convexity parameter
        NT kappa; // Condition number

    parameters() :
        alpha(NT(1)),
        order(2),
        L(NT(1)),
        m(NT(1)),
        kappa(1)
    {};

    parameters(
        NT alpha_,
        unsigned int order_) :
        alpha(alpha_),
        order(order_),
        L(alpha_),
        m(alpha_),
        kappa(1)
    {};
  };


    template
    <
        typename Point
    >
    struct GradientFunctor {
        typedef typename Point::FT NT;
        typedef std::vector<Point> pts;

        parameters<NT> &params;

        GradientFunctor(parameters<NT> &params_) : params(params_) {};

        // The index i represents the state vector index
        Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
            if (i == params.order - 1) {
                return (-params.alpha) * xs[0]; // returns - a*x
            } else {
                return xs[i + 1]; // returns derivative
        }
    }

    Point operator()(Point const &x){
      Point y = (-params.alpha) * x;
      return y;
    }
  };


    template
    <
        typename Point
    >
    struct FunctionFunctor {
        typedef typename Point::FT NT;

        parameters<NT> &params;

        FunctionFunctor(parameters<NT> &params_) : params(params_) {};

        NT operator() (Point const& x) const {
            return 0.5 * params.alpha * x.dot(x);
        }
  };

};

struct IsotropicLinearFunctor {

    // Exponential Density
    template <
        typename NT
    >
    struct parameters {
        NT alpha;
        unsigned int order;
        NT L; // Lipschitz constant of gradient
        NT m; // Strong-convexity constant
        NT kappa; // Condition number

    parameters() :
        alpha(NT(1)),
        order(1),
        L(0),
        m(0),
        kappa(1)
    {};

    parameters(NT alpha_, unsigned int order_) :
        alpha(alpha_),
        order(order),
        L(0),
        m(0),
        kappa(1)
    {}
  };

    template
    <
        typename Point
    >
    struct GradientFunctor {
        typedef typename Point::FT NT;
        typedef std::vector<Point> pts;

        parameters<NT> &params;

        GradientFunctor(parameters<NT> &params_) : params(params_) {};

        // The index i represents the state vector index
        Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
            if (i == params.order - 1) {
                Point y = Point::all_ones(xs[0].dimension());
                y = (- params.alpha) * y;
                return y;
            } else {
                return xs[i + 1]; // returns derivative
        }
    }

  };

    template
    <
        typename Point
    >
    struct FunctionFunctor {
        typedef typename Point::FT NT;

        parameters<NT> &params;

        FunctionFunctor(parameters<NT> &params_) : params(params_) {};

        NT operator() (Point const& x) const {
            return params.alpha * x.sum();
        }

  };

};


struct ExponentialFunctor {

  // Sample from linear program c^T x (exponential density)
  template <
      typename NT,
      typename Point
  >
  struct parameters {
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number
    Point c; // Coefficients of LP objective
    NT a; // Inverse variance

    parameters(Point c_) : order(2), L(1), m(1), kappa(1), c(c_), a(1.0) {};
    parameters(Point c_, NT a_) : order(2), L(1), m(1), kappa(1), c(c_), a(a_) {};

  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y(params.c);
        return (-params.a) * y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }

  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      return params.a * x.dot(params.c);
    }

  };

};


struct GaussianFunctor {

  template <
      typename NT,
      typename Point
  >
  struct parameters {
    Point x0;
    NT a;
    NT eta;
    unsigned int order;
    NT L; // Lipschitz constant for gradient
    NT m; // Strong convexity constant
    NT kappa; // Condition number

    parameters(Point x0_, NT a_, NT eta_) :
        x0(x0_), a(a_), eta(eta_), order(2), L(2 * a_), m(2 * a_), kappa(1) {};

  };

  template
  <
      typename Point
  >
  struct GradientFunctor {
    typedef typename Point::FT NT;
    typedef std::vector<Point> pts;

    parameters<NT, Point> &params;

    GradientFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    Point operator() (unsigned int const& i, pts const& xs, NT const& t) const {
      if (i == params.order - 1) {
        Point y = (-2.0 * params.a) * (xs[0] - params.x0);
        return y;
      } else {
        return xs[i + 1]; // returns derivative
      }
    }
    Point operator()(Point const&x){
      Point y = (-2.0 * params.a) * (x - params.x0);
      return y;
    }
  };

  template
  <
    typename Point
  >
  struct FunctionFunctor {
    typedef typename Point::FT NT;

    parameters<NT, Point> &params;

    FunctionFunctor(parameters<NT, Point> &params_) : params(params_) {};

    // The index i represents the state vector index
    NT operator() (Point const& x) const {
      Point y = x - params.x0;
      return params.a * y.dot(y);
    }

  };

  template
<
  typename Point
>
struct HessianFunctor {
  typedef typename Point::FT NT;

  parameters<NT, Point> &params;

  HessianFunctor(parameters<NT, Point> &params_) : params(params_) {};

  // The index i represents the state vector index
  Point operator() (Point const& x) const {
    return (2.0 * params.a) * Point::all_ones(x.dimension());
  }

};

};

#endif
