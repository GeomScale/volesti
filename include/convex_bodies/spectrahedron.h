//
// Created by panagiotis on 25/6/2019.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "Eigen"
#include <vector>
#include <Eigen/Eigen>
#include <limits>


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;


/**
 * A linear matrix inequality A_0 + sum(x * F_i)
 */
class LMI {
    MT A0;

    // the matrices A_i, i>0
    std::vector<MT> matrices;

    typedef std::vector<MT>::iterator Iter;

public:
    LMI() {};

    LMI(std::vector<MT>& matrices) {
        this->A0 = matrices[0];

        for (int i=1 ; i<matrices.size() ; i++)
            this->matrices.push_back(matrices[i]);
    }

    /**
     * Evaluate the lmi for vector x
     *
     * @param x
     * @return
     */
    MT evaluate(VT& x) {
       MT res = A0;
       int i = 0;

       for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++)
           res += x(i) * (*iter);
    }

    /**
     * Evaluate the lmi for vector x without taking int account matrix A0
     *
     * @param x
     * @return
     */
    MT evaluateWithoutA0(VT& x) {
        long dim = A0.rows();
        MT res;
        res.setZero(dim, dim);
        int i = 0;

        for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++)
            res += x(i) * (*iter);
    }

    const MT& getA0() const {
        return A0;
    }

};


class Spectrahedron {
    /**
     * The collection of matrices that constitute the linear matrix
     * inequality describing the spectrahedron
     */
    LMI lmi;

    double maxDouble = std::numeric_limits<double>::max();
    double minDouble = std::numeric_limits<double>::lowest();

public:

    Spectrahedron() {};

    Spectrahedron(LMI& lmi) {
        this->lmi = lmi;
    };

    const LMI& getLMI() const {
        return lmi;
    }

    /**
     * Compute the intersection of a 1D line and the spectrahedron by finding the
     * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
     *
     * @param position
     * @param direction
     * @return (minimum positive eigenvalue, maximum negative eigenvalue)
     */
    std::pair<double, double> boundaryOracle(VT& position, VT& direction) {
        MT A = lmi.evaluate(position);

        MT B = -lmi.evaluateWithoutA0(direction);

        Eigen::GeneralizedEigenSolver<MT> ges;
        ges.compute(A, B, false);

        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eigenvalues = ges.eigenvalues(); //TODO not like this check Eigen may divide by zero

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;

        for (int i=0 ; i<eigenvalues.rows() ; i++) {
            if (eigenvalues(i).imag() != 0)
                continue;

            if (eigenvalues(i).real() > 0 && eigenvalues(i).real() < lambdaMinPositive)
                lambdaMinPositive = eigenvalues(i).real();
            if (eigenvalues(i).real() < 0 && eigenvalues(i).real() > lambdaMaxNegative)
                lambdaMaxNegative = eigenvalues(i).real();
        }

        return {lambdaMinPositive, lambdaMaxNegative};
    }

    /**
    * Compute the intersection of a 1D line and the spectrahedron by finding the
    * generalized eigenvalues of (LMI(x), A0 - LMI(direction))
    *
    * Take also in account the halfspace ax<=b
     *
    * @param position
    * @param direction
    * @return (minimum positive eigenvalue, maximum negative eigenvalue)
    */
    std::pair<double, double> boundaryOracle(VT& position, VT& direction, VT& a, double b) {
        MT A = lmi.evaluate(position);

        MT B = -lmi.evaluateWithoutA0(direction);

        Eigen::GeneralizedEigenSolver<MT> ges;
        ges.compute(A, B, false);

        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eigenvalues = ges.eigenvalues(); //TODO not like this check Eigen may divide by zero

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;

        for (int i=0 ; i<eigenvalues.rows() ; i++) {
            if (eigenvalues(i).imag() != 0)
                continue;

            if (eigenvalues(i).real() > 0 && eigenvalues(i).real() < lambdaMinPositive)
                lambdaMinPositive = eigenvalues(i).real();
            if (eigenvalues(i).real() < 0 && eigenvalues(i).real() > lambdaMaxNegative)
                lambdaMaxNegative = eigenvalues(i).real();
        }

        double lambda = (b - a.dot(position)) / a.dot(direction);

        if (lambda > 0 && lambda < lambdaMinPositive)
            lambdaMinPositive = lambda;
        if (lambda < 0 && lambda > lambdaMaxNegative)
            lambdaMaxNegative = lambda;

        return {lambdaMinPositive, lambdaMaxNegative};
    }

};

#endif //VOLESTI_SPECTRAHEDRON_H
