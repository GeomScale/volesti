// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ODE_SOLVERS_BASIS_HPP
#define ODE_SOLVERS_BASIS_HPP

enum BasisType {
  DERIVATIVE = 0,
  FUNCTION = 1,
  INTEGRAL = 2
};


template <typename NT>
struct PolynomialBasis {
  BasisType basis_type;

  PolynomialBasis(BasisType basis_type_) : basis_type(basis_type_) {}

  NT operator() (NT t, const NT t0, unsigned int j, const unsigned int ord) const {
    switch (basis_type) {
      case FUNCTION:
        return pow(t - t0, NT(j));
      case DERIVATIVE:
        return NT(j) * pow(t - t0, NT(j - 1));
      case INTEGRAL:
        return pow(t - t0, NT(j + 1)) / NT(j + 1);
    }
  }

};

template <typename NT>
struct Polynomial {

  BasisType basis_type;
  PolynomialBasis<NT> basis;
  std::vector<NT> coeffs;
  NT result;
  unsigned int ord;

  Polynomial(std::vector<NT> coeffs_, BasisType basis_type_) :
    basis_type(basis_type_), coeffs(coeffs_), basis(basis_type_) {
      ord = coeffs.size();
    }

  NT operator() (NT t, NT t0, unsigned int j=-1, unsigned int ord=-1) {
    result = NT(0);

    for (unsigned int i = 0; i < ord; i++) {
        result += coeffs[i] * basis(t, t0, i, ord);
    }

    return result;
  }

  static std::vector<NT> convolve(std::vector<NT> const& p, std::vector<NT> const& q) {
    // Performs direct convolution of p and q (less round-off error than FFT)
    unsigned int n = p.size();
    unsigned int m = q.size();
    unsigned int r = n + m - 1;

    std::vector<NT> result(r, NT(0));

    unsigned int j, k;

    for (unsigned int i = 0; i < r; i++) {
      j = (i >= m - 1)? i - (m - 1) : 0;
      k = (i <  n - 1)? i : n - 1;
      for (unsigned int z = j; z <= k; z++) result[i] += (p[z] * q[i - z]);
    }

    return result;
  }

  static std::vector<NT> multi_convolve(std::vector<std::vector<NT>> &seq) {

    std::vector<NT> result;
    result = seq[0];

    for (unsigned int i = 1; i < seq.size(); i++) {
      result = convolve(result, seq[i]);
    }

    return result;

  }


};

template <typename NT, typename bfunc>
struct RationalFunctionBasis {
  bfunc p, q;
  bfunc grad_p, grad_q;
  NT reg = 1e-6;
  BasisType basis_type;
  NT num, den, grad_num, grad_den;

  RationalFunctionBasis(bfunc num, bfunc grad_num, bfunc den, bfunc grad_den, BasisType basis_type_) :
    p(num), grad_p(grad_num), q(den), grad_q(grad_den), basis_type(basis_type_) {};

  NT operator()(NT t, NT t0, unsigned int j, unsigned int ord) {

    switch (basis_type) {
      case FUNCTION:
        num = p(t, t0, j, ord);
        den = q(t, t0, j, ord);
        if (std::abs(den) < reg) den += reg;
        return num / den;
      case DERIVATIVE:
        num = p(t, t0, j, ord);
        grad_num = grad_p(t, t0, j, ord);
        den = q(t, t0, j, ord);
        grad_den = grad_q(t, t0, j, ord);
        if (std::abs(den * den) < reg) den += reg;
        return (grad_num  / den)  - (grad_den * num) / den;
      case INTEGRAL:
        throw true;
    }
  }

};

template <typename NT, typename VT>
struct LagrangePolynomial {

  VT coeffs;
  VT nodes;
  int basis = -1;

  LagrangePolynomial() {}

  int order() const {
    return nodes.rows();
  }

  void set_nodes(VT nodes_) {
    nodes = nodes_;
  }

  void set_basis(int basis_) {
    basis = basis_;
  }

  void set_coeffs(VT coeffs_) {
    coeffs = coeffs_;
  }

  NT operator() (NT t) const {

    int j = (int) round((order() * acos(t)) / M_PI - 0.5);
    NT temp = cos((j+0.5) * M_PI / order());

    if (abs(temp - nodes(j)) < 1e-6) {
      if (basis == -1) return coeffs(j);
      else return NT(1);
    } else {
      throw true;
      return NT(-1);
    }

  }
};

template <typename T>
void degree_doubling_chebyshev(std::vector<T> &coeffs,
  std::vector<T> &result) {
  unsigned int N = coeffs.size() - 1;

  for (int i = 0; i <= 2 * N; i++) {
    if (i > N) result[i] = coeffs[i - N];
    else if (i == N) result[i] = 2 * coeffs[0];
    else result[i] = coeffs[N - i];
  }

}

#endif
