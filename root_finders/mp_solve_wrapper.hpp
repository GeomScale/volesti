// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MP_SOLVE_WRAPPER_HPP
#define MP_SOLVE_WRAPPER_HPP

template <typename NT>
std::vector<std::pair<NT, NT>> mpsolve(std::vector<NT> &coeffs, bool positive_real=false) {

  long n = (long) coeffs.size();

  while (std::abs(coeffs[n-1]) < NT(1e-9)) {
    n--;

  }


  mps_monomial_poly *p;
  mps_context *s;

  s = mps_context_new ();
  p = mps_monomial_poly_new (s, n-1);

  mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);

  for (long i = 0; i < n; i++) {
      mps_monomial_poly_set_coefficient_d (s, p, i, coeffs[i], 0);
  }


  /* Set the input polynomial */
  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

  /* Allocate space to hold the results. We check only floating point results
  * in here */
  cplx_t *results = cplx_valloc (n-1);

  /* Actually solve the polynomial */
  mps_mpsolve (s);

  /* Save roots computed in the vector results */
  mps_context_get_roots_d (s, &results, NULL);

  std::vector<std::pair<NT, NT>> results_vector;

  NT real, im;

  for (long i = 0; i < n - 1; i++) {
    real = (NT) cplx_Re(*results);
    im = (NT) cplx_Im(*results);

    #ifdef VOLESTI_DEBUG
      std::cout << real << " + " << im << "i" << std::endl;
    #endif
    results++;
    if (positive_real) {
      if (real > 0 && std::abs(im) < 1e-8) {
        results_vector.push_back(std::make_pair(real, im));
      }
    } else {
      results_vector.push_back(std::make_pair(real, im));
    }
  }

  return results_vector;
}

#endif
