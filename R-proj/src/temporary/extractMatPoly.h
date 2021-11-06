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


#ifndef EXTRACTMATPOLY_H
#define EXTRACTMATPOLY_H

// Take a H or a V-polytope and return a numerical matrix in ine or ext format respectively
template <class Polytope>
Rcpp::NumericMatrix extractMatPoly(Polytope P) {

    typedef typename Polytope::MT 	MT;

    MT Mat(P.get_mat().rows(), P.dimension()+1);
    Mat << P.get_vec(), P.get_mat();

    return Rcpp::wrap(Mat);
}

#endif
