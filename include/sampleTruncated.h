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

#ifndef SAMPLETRUNCATED_H
#define SAMPLETRUNCATED_H

#include "TruncatedNormal2.h"

using namespace MCQMCIntegration;

void mvrandn(int N) {

    const uint32_t s = 5;
    double LOWER[5] = {-INFINITY, -INFINITY, -INFINITY, -INFINITY, -INFINITY};
    double UPPER[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double SIGMA[5][5] = {
            {1.666667e+00, 8.333333e-01, 8.333333e-01, 8.333333e-01, 8.333333e-01},
            {8.333333e-01, 1.666667e+00, 8.333333e-01, 8.333333e-01, 8.333333e-01},
            {8.333333e-01, 8.333333e-01, 1.666667e+00, 8.333333e-01, 8.333333e-01},
            {8.333333e-01, 8.333333e-01, 8.333333e-01, 1.666667e+00, 8.333333e-01},
            {8.333333e-01, 8.333333e-01, 8.333333e-01, 8.333333e-01, 1.666667e+00}
    };

    vector<double> lower(LOWER, LOWER+s);
    vector<double> upper(UPPER, UPPER+s);
    vector< vector<double> > sigma(s);
    for (uint32_t i = 0; i < s; i++) {
        for (uint32_t j = 0; j < s; j++) {
            sigma[i].push_back(SIGMA[i][j]);
        }
    }
    // Monte-Carlo Method
    const TruncatedNormalResult result
            = truncatedNormalMC(lower,
                                upper,
                                sigma,
                                10000, // number of trials.
                                1000, // number of samples per a trial.
                                0.99);
    if (!result.success) {
        cout << "calculation failed" << endl;
        return;
    }
    // show result
    cout << "probability:" << result.probability << endl;
    cout << "absoluteError:" << result.absoluteError << endl;
    cout << "relativeError:" << result.relativeError << endl;
    cout << "theoreticalUpperBound:" << result.upperBound << endl;
    return;

}


#endif