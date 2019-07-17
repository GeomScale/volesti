//
// Created by panagiotis on 9/7/2019.
//

#ifndef VOLESTI_SDP_GENERATOR_H
#define VOLESTI_SDP_GENERATOR_H

#include "spectrahedron.h"
#include "Eigen"
#include "sdp_problem.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;

// m must be even
template <class Point>
optimization::sdp_problem<Point> generateSDP(int n, int m) {
    VT obj(n);
    srand (time(NULL));
    
    if (m % 2 != 0) throw 1;

    for (int i = 0; i < n; i++) {
        obj(i) = -5 + ((double)rand() / RAND_MAX)*30;
    }

    MT ones(m, m);
    for (int i=0 ; i<ones.rows() ; i++)
        for (int j=0 ; j<ones.cols() ; j++)
            ones(i, j) = 1;

    MT M = 2* Eigen::MatrixXd::Random(m,m) - ones;
    MT I = Eigen::MatrixXd::Identity(m, m);

    std::vector<MT> matrices(n + 1);
    matrices[0] = -(M * M.transpose()) - I;

    ones.resize(m/2, m/2);
    for (int i=0 ; i<ones.rows() ; i++)
        for (int j=0 ; j<ones.cols() ; j++)
            ones(i, j) = 1;

    for (int i=1 ; i<=n ; i++) {
        M =  2* Eigen::MatrixXd::Random(m/2,m/2) - ones;
//        MT F = M + M.transpose();
        MT F(m/2, m/2);

        for (int j=0 ; j<m/2 ; j++)
            for (int k=j ; k<m/2 ; k++)
                F(j,k) = M(j,k);

        MT FF = F + F.transpose();

        for (int j=0 ; j<m/2 ;j++)
            FF(j,j) -= M(j,j);

        MT A;
        A.setZero(m, m);
        A.topLeftCorner(m/2, m/2) = FF;
        A.bottomRightCorner(m/2, m/2) = -FF;
        matrices[i] = A;
    }

    LMI lmi(matrices);
    Spectrahedron spectrahedron(lmi);

    return optimization::sdp_problem<Point>(spectrahedron, obj);
}

#endif //VOLESTI_SDP_GENERATOR_H
