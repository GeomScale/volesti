/*
GeomScale Project

Copyright (c) 2020
  Vissarion Fisikopoulos
  Apostolos Chalkis
  Elias Tsigaridas
  Marios Papachristou

Contributed and/or modified by Marios Papachristou,
as part of Google Summer of Code 2020 program.

VolEsti is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

VolEsti is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

See the file COPYING.LESSER for the text of the GNU Lesser General
Public License.  If you did not receive this file along with HeaDDaCHe,
see <http://www.gnu.org/licenses/>.
*/

#ifndef MP_SOLVE_WRAPPER_H
#define MP_SOLVE_WRAPPER_H

template <typename NT>
std::vector<std::pair<NT, NT>> mpsolve(NT t0, std::vector<NT> &coeffs, bool positive_real=false) {

  long n = (long) coeffs.size();
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
