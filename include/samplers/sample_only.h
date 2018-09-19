// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.


#ifndef SAMPLE_ONLY_H
#define SAMPLE_ONLY_H

template <class Point, typename NT, class PointList, class Polytope, class UParameters, class GParameters>
void sampling_only(PointList &randPoints, Polytope &P, unsigned int walk_len,
                   unsigned int rnum, bool gaussian, NT a, Point internal_point,
                   UParameters var1, GParameters var2) {

    typedef typename UParameters::RNGType RNGType;
    unsigned int n = var1.n;
    Point p = internal_point;
    Point q = get_point_on_Dsphere<RNGType, Point>(n, var1.che_rad);
    p=p+q;
    rand_point_generator(P, p, 1, 50 * n, randPoints, var1);

    randPoints.clear();
    if (!gaussian){
        rand_point_generator(P, p, rnum, walk_len, randPoints, var1);
    } else {
        rand_gaussian_point_generator(P, p, rnum, walk_len, randPoints, a, var2);
    }

}

#endif
