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
template <class T>
Rcpp::NumericMatrix extractMatPoly(T &P) {


    Eigen::MatrixXd A = P.get_mat();
    Eigen::VectorXd b = P.get_vec();
    int n = P.dimension(), m = A.rows(), i, j;

    Rcpp::NumericMatrix Mat(m+1, n+1);

    Mat(0,0) = m; Mat(0,1) = n+1;
    for (i=2; i<n+1; i++){
        Mat(0,i) = 0;
    }

    for (i=1; i<m+1; i++){
        Mat(i,0) = b(i-1);
        for (j=1; j<n+1; j++){
            Mat(i,j) = A(i-1,j-1);
        }
    }

    return Mat;
}

#endif
