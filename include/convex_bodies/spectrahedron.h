//
// Created by panagiotis on 25/6/2019.
//

#ifndef VOLESTI_SPECTRAHEDRON_H
#define VOLESTI_SPECTRAHEDRON_H

#include "Eigen"
#include <vector>
#include <Eigen/Eigen>
#include <limits>

const double ZERO = 0.0000000001;

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

    LMI(const LMI& lmi) {
        this->A0 = lmi.A0;
        this->matrices = lmi.matrices;
    }

    /**
     * Evaluate the lmi for vector x
     *
     * @param x
     * @return
     */
    MT evaluate(const VT& x) {
       MT res = A0;
       int i = 0;

       for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++)
           res += x(i) * (*iter);

       return res;
    }

    bool isSingular(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() == 0)
                return true;

        return false;
    }

    bool isSingular(VT& x, double approx) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (abs(eivals(i).real()) <= abs(approx))
                return true;

        return false;
    }

    bool isNegativeDefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() >= 0)
                return false;

        return true;
    }

    bool isNegativeSemidefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() > 0)
                return false;

        return true;
    }

    bool isPositiveSemidefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() < 0)
                return false;

        return true;
    }

    bool isPositiveDefinite(VT& x) {
        MT mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

        for (int i = 0; i < eivals.rows(); i++)
            if (eivals(i).real() <= 0)
                return false;

        return true;
    }


    bool isPositiveSemidefinite(VT& x, MT& mt, VT& minEigenVector, double& minEigenvalue) {
        bool res = true;

        mt = evaluate(x);

        Eigen::EigenSolver<MT> solver;
        solver.compute(mt);
        Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
        minEigenVector.resize(solver.eigenvectors().col(0).rows());

        for (int i=0 ; i<minEigenVector.rows() ; i++) {
            minEigenVector(i) = solver.eigenvectors().col(0)(i).real();

        }

        double min = eivals(0).real() / minEigenVector.norm();
        int minIndex = 0;

        for (int i = 0; i < eivals.rows(); i++) {
            for (int j=0 ; j<minEigenVector.rows() ; j++) {
                minEigenVector(j) = solver.eigenvectors().col(i)(j).real();
            }

            if (eivals(i).real() < 0)
                res = false;

            if (eivals(i).real() / minEigenVector.norm() < min) {
                min = eivals(i).real() / minEigenVector.norm();
                minIndex = i;
            }
        }


        minEigenvalue = eivals(minIndex).real()  / minEigenVector.norm() ;

        for (int i=0 ; i<minEigenVector.rows() ; i++) {
            minEigenVector(i) = solver.eigenvectors().col(minIndex)(i).real();

        }
        minEigenVector.normalize();

        return res;
    }


    int getDim() {
        return matrices.size();
    }

    const std::vector<MT>& getMatrices() const {
        return matrices;
    }

    /**
     * Evaluate the lmi for vector x without taking int account matrix A0
     *
     * @param x
     * @return
     */
    MT evaluateWithoutA0(const VT& x) {
        long dim = A0.rows();
        MT res;
        res.setZero(dim, dim);
        int i = 0;

        for (Iter iter=matrices.begin() ; iter!=matrices.end() ; iter++, i++)
            res += x(i) * (*iter);

        return res;
    }

    const MT& getAi(int i) const {
        return matrices[i];
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

    void changeInequalityDirection() {
        A0 = -1 * A0;

        for (int i=0 ; i<matrices.size() ; i++)
            matrices[i] = -1 * matrices[i];
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

    Spectrahedron(const Spectrahedron& spectrahedron) {
        LMI lmi;
        this->lmi = LMI(spectrahedron.getLMI());
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
    std::pair<double, double> boundaryOracle(const VT& position, const VT& direction) {
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
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
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
    std::pair<double, double> boundaryOracle(const VT& position, const VT& direction, const VT& a, double b) {
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
//        if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//        if (lambdaMaxNegative > -ZERO) lambdaMaxNegative = 0;
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

    void changeInequalityDirection() {
        lmi.changeInequalityDirection();
    }

    bool isSingular(VT& x) {
        return lmi.isSingular(x);
    }


    bool isSingular(VT& x, double approx) {
        return lmi.isSingular(x, approx);
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
