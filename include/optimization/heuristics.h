//
// Created by panagiotis on 21/7/2019.
//

#ifndef VOLESTI_HEURISTICS_H
#define VOLESTI_HEURISTICS_H

#include "Eigen"

namespace optimization {

    /**
     * For the Eigen library
     */
    typedef double NT_MATRIX;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, 1> VT;


    template <class Point>
    MT sampledCovarianceMatrix(std::list<Point> &points) {
        int dim = points.front().dimension();
        double pointsNum = points.size();

        VT y;
        y.setZero(dim);

        for (auto p : points)
            y += p.getCoefficients() / pointsNum;

        // compute Y
        MT Y;
        Y.setZero(dim, dim);
        VT temp(dim);

        for (auto pit = points.begin(); pit != points.end(); pit++) {
            temp = pit->getCoefficients() - y;
            Y = Y + ((temp * temp.transpose()) / (pointsNum - 1));
        }

        try {
            Y = Y.sqrt();
        }
        catch (int e) {
            throw e;
        }

        return Y;
    }
}

#endif //VOLESTI_HEURISTICS_H
