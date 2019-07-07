//
// Created by panagiotis on 25/6/2019.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "Eigen"
#include <vector>
#include <Eigen/Eigen>
#include <limits>

const double ZERO = 0.000000000001;

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

    LMI(LMI& lmi) {
        this->A0 = lmi.A0;
        this->matrices = lmi.matrices;
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

       return res;
    }

    bool isNegativeDefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
        std::cout << eivals << "\n" <<  solver.eigenvectors() << "\n";

        for (int i = 1; i < eivals.rows(); i++)
            if (eivals(i).real() > 0)
                return false;

        return true;
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

        return res;
    }

    const MT& getA0() const {
        return A0;
    }

    void setA0(MT& A0) {
        this->A0 = A0;
    }

    void addMatrix(MT& matrix) {
        matrices.push_back(matrix);
    }

    void print() {
        std::cout << "F0" << "\n" << A0 << "\n";
        int i = 1;

        for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++) {
            std::cout << "F" << i << "\n";
            std::cout << *iter << "\n";
        }
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

    Spectrahedron() {}

    Spectrahedron(Spectrahedron& spectrahedron) {
        this->lmi = LMI(spectrahedron.lmi);
    }

    Spectrahedron(LMI& lmi) {
        this->lmi = lmi;
    }

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


        Eigen::GeneralizedEigenSolver<MT> ges(A,B);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();


        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;


        for (int i=0 ; i<alphas.rows() ; i++) {
            if (betas(i) == 0) //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative =lambda;
        }

        // for numerical stability
        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive ==  maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
        if (lambdaMaxNegative == minDouble) lambdaMaxNegative = 0;


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

        Eigen::GeneralizedEigenSolver<MT> ges(A, B);

        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
        VT betas = ges.betas();

        double lambdaMaxNegative = minDouble;
        double lambdaMinPositive = maxDouble;

        for (int i=0 ; i<alphas.rows() ; i++) {

            if (betas(i) == 0)  //TODO WARNING do what here?
                continue;

            double lambda = alphas(i).real() / betas(i);

            if (lambda > 0 && lambda < lambdaMinPositive)
                lambdaMinPositive = lambda;
            if (lambda < 0 && lambda > lambdaMaxNegative)
                lambdaMaxNegative =lambda;
        }


        // for numerical stability
        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
        if (lambdaMinPositive ==  maxDouble) lambdaMinPositive = 0; //TODO b must be too small..
        if (lambdaMaxNegative == minDouble) lambdaMaxNegative = 0;


        // check the cutting plane
        double lambda = (b - a.dot(position)) / a.dot(direction);
        if (lambda > 0 && lambda < lambdaMinPositive)
            lambdaMinPositive = lambda;
        if (lambda < 0 && lambda > lambdaMaxNegative)
            lambdaMaxNegative = lambda;


        return {lambdaMinPositive, lambdaMaxNegative};
    }

    void print() {
        this->lmi.print();
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
