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

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

template <typename NT, class bfunc>
NT newton_raphson(NT t0, bfunc f, bfunc grad_f, NT rtol) {
  NT t, t_prev, err;
  t = t0;

  do {
    t = t_prev - (f(t_prev) / grad_f(t_prev));
    if (t_prev != 0) {
      err = std::abs(t - t_prev) / t_prev;
    } else {
      err = std::abs(t - t_prev);
    }

    t_prev = t;
  } while (err > rtol);

  return t;
}

#endif
