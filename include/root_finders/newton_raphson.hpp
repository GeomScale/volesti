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

#ifndef NEWTON_RAPHSON_HPP
#define NEWTON_RAPHSON_HPP

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
