//
// Created by panagiotis on 7/7/2019.
//

#ifndef VOLESTI_LMI_STRICT_FEASIBILITY_H
#define VOLESTI_LMI_STRICT_FEASIBILITY_H

#include "Eigen"
#include "spectrahedron.h"
#include <limits>

//typedef double NT_MATRIX;
//typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, Eigen::Dynamic> MT;
//typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, 1> VT;
//
//
//double boundaryOracle(LMI& lmi, MT& A, VT& position, VT& direction) {
//    Eigen::EigenSolver<MT> solver;
//    solver.compute(A);
//    Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();
//    bool ApositiveDefinite = true, AnegativeDefinite = true;
//
//    for (int i = 0; i < eivals.rows(); i++) {
//        if (eivals(i).real() > 0) AnegativeDefinite = false;
//        if (eivals(i).real() < 0) ApositiveDefinite = false;
//    }
//
//    MT B = -lmi.evaluateWithoutA0(direction);
//
//    solver.compute(A);
//    eivals = solver.eigenvalues();
//    bool BpositiveDefinite = true, BnegativeDefinite = true;
//
//    for (int i = 0; i < eivals.rows(); i++) {
//        if (eivals(i).real() > 0) BnegativeDefinite = false;
//        if (eivals(i).real() < 0) BpositiveDefinite = false;
//    }
//
//    std::cout << ApositiveDefinite << " " << AnegativeDefinite << " " << BpositiveDefinite << " " << BnegativeDefinite << "\n";
//    Eigen::GeneralizedEigenSolver<MT> ges(A,B);
//    Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType alphas = ges.alphas();
//    VT betas = ges.betas();
//
//    double lambdaMinPositive = std::numeric_limits<double>::max();
//    double lambdaMaxNegative = std::numeric_limits<double>::lowest();
//
//    std::cout<<ges.eigenvalues() <<"\n";
//
//    for (int i=0 ; i<alphas.rows() ; i++) {
//        if (betas(i) == 0) //TODO WARNING do what here?
//            continue;
//
//        double lambda = alphas(i).real() / betas(i);
//
//        if (lambda > 0 && lambda < lambdaMinPositive)
//            lambdaMinPositive = lambda;
//        if (lambda < 0 && lambda > lambdaMaxNegative)
//            lambdaMaxNegative = lambda;
//    }
//
//    // for numerical stability
//    if (lambdaMinPositive < ZERO) lambdaMinPositive = 0;
//    if (lambdaMinPositive ==  std::numeric_limits<double>::max()) lambdaMinPositive = 0; //TODO b must be too small..
//
//    return -lambdaMaxNegative;
//}
//
//void getFeasibilityVector(LMI& lmi, MT& A, VT& minEigenVector, double minEigenValue, VT& feasibilityVector, double& feasibilityDistance) {
//    int dim = lmi.getDim();
//    feasibilityVector.resize(dim);
//    minEigenVector.normalize();
//
//    VT gradient(dim);
//    for (int i=0 ; i<dim ; i++) {
//        MT Ai = lmi.getAi(i);
//        VT temp = minEigenVector.transpose() * Ai;
//        gradient(i) = temp.dot(minEigenVector);
//    }
//
//    double gradientNormSquared = 0;
//    for (int i=0 ; i<dim ; i++) {
//        gradientNormSquared += gradient(i) * gradient(i);
//    }
//
//    feasibilityVector = -(minEigenValue * gradient) / gradientNormSquared;
//    feasibilityDistance = std::abs(minEigenValue / std::sqrt(gradientNormSquared));
//}
//
//bool DBmax(LMI& lmi, VT& point, int iterationsNum, double feasibilityDistanceTolerance, double movementTolerance) {
//
//    bool success;
//    VT s_plus, s_minus, n_plus, n_minus;
//    int dim = point.rows();
//
//    for (int iter=0 ; iter<iterationsNum ; iter++) {
//        VT minEigenVector;
//        double minEigenValue;
//        MT evalLMI;
//        s_plus.setZero(dim);
//        s_minus.setZero(dim);
//        n_plus.setZero(dim);
//        n_minus.setZero(dim);
//        success = true;
//
//        if (!lmi.isPositiveSemidefinite(point, evalLMI, minEigenVector, minEigenValue)) {
//            success = false;
//            VT feasibilityVector;
//            double feasibilityDistance;
//            getFeasibilityVector(lmi, evalLMI, minEigenVector, minEigenValue, feasibilityVector,
//                                 feasibilityDistance);
//            std::cout << iter << " " << feasibilityDistance << "\n";
//            std::cout << feasibilityVector.transpose() << "\n" << point.transpose() << "  wddddddddddd\n";
//
//            if (feasibilityDistance > feasibilityDistanceTolerance) {
//
//                for (int i = 0; i < feasibilityVector.rows(); i++) {
//                    if (feasibilityVector(i) > 0) {
//                        n_plus(i) = n_plus(i) + 1;
//                        if (feasibilityVector(i) > s_plus(i)) s_plus(i) = feasibilityVector(i);
//                    } else if (feasibilityVector(i) < 0) {
//                        n_minus(i) = n_minus(i) + 1;
//                        if (feasibilityVector(i) < s_minus(i)) s_minus(i) = feasibilityVector(i);
//                    }
//                }
//            }
//        }
//
//        if (success) return true;
//
//        VT t(dim);
//
//        for (int i = 0; i < point.rows(); i++) {
//            if (n_plus(i) == n_minus(i))
//                t(i) = (s_plus(i) + s_minus(i)) / 2;
//            else if (n_plus(i) > n_minus(i))
//                t(i) = s_plus(i);
//            else
//                t(i) = s_minus(i);
//        }
//
//        double tNorm = 0;
//        for (int i=0 ; i<dim ; i++) {
//            tNorm += t(i) * t(i);
//        }
//        tNorm = std::sqrt(tNorm);
//
//        if (tNorm < movementTolerance)
//            return false;
//
//        point = point + t;
//    } /* for (int iter=0 ; iter<iterationsNum ; iter++) { */
//assert(false);
//    return  false;
//}
//
//bool basicConsensusMethod(LMI& lmi, VT& x, int iterationsNum) {
//    bool infeasible = true;
//    int iter = 0;
//    VT s;
//    int dim = x.rows();
//
//    while (iter++ < iterationsNum && infeasible) {
//        s.setZero(dim);
//        VT minEigenVector;
//        double minEigenValue;
//        MT evalLMI;
//        VT feasibilityVector;
//        double feasibilityDistance;
//
//        if (!lmi.isPositiveSemidefinite(x, evalLMI, minEigenVector, minEigenValue)) {
//            getFeasibilityVector(lmi, evalLMI, minEigenVector, minEigenValue, feasibilityVector,
//                                 feasibilityDistance);
//
//            double l_plus = boundaryOracle(lmi, evalLMI, x, feasibilityVector);
//            std::cout << l_plus <<"\n";
////            VT z = x + 0.5 * feasibilityVector;
//            x = x + 0.5 * (2*l_plus )* feasibilityVector;
//        }
//        else
//            return true;
//
//    } /* while (iter++ > iterationsNum && infeasible) { */
//
//    return !infeasible;
//}
//
//VT strict_feasible_point(Spectrahedron& spectrahedron, double distance_tolerance, double movement_tolerance, int iterationsNum) {
//    LMI lmi;
//    lmi = spectrahedron.getLMI();
//    int dim = lmi.getDim();
//    VT x(dim);
//
//    DBmax(lmi, x, iterationsNum, distance_tolerance, movement_tolerance);
//    basicConsensusMethod(lmi, x, iterationsNum);
//    return x;
//}


#endif //VOLESTI_LMI_STRICT_FEASIBILITY_H
