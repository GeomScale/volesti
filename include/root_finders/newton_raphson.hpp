// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NEWTON_RAPHSON_HPP
#define NEWTON_RAPHSON_HPP

/// @brief Function implementing the Newton-Raphson numerical method
/// @tparam NT Number type
/// @tparam func Function type
template <typename NT, class func>
std::pair<NT,bool> newton_raphson(NT t0, func f, func grad_f, const NT rtol,
  const NT reg=0, const unsigned int max_tries=1000000) {
  NT t, t_prev, err;
  NT y, y_prime;
  t = t0;
  unsigned int tries = 0;

  do {
    tries++;
    y = f(t_prev);
    y_prime = grad_f(t_prev);

    if (std::abs(y_prime) < rtol) y_prime += reg;

    t = t_prev - y / y_prime;
    if (t_prev != 0) {
      err = std::abs(t - t_prev) / t_prev;
    } else {
      err = std::abs(t - t_prev);
    }

    t_prev = t;

    if (tries > max_tries) break;

  } while (err > rtol);

  return std::make_pair(t, tries > max_tries);

}

#endif
