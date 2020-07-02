#ifndef SOS_UTILS_H
#define SOS_UTILS_H

#include <Eigen/Dense>
#include <iostream>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cxxtimer.hpp>
#include "spdlog/spdlog.h"


typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> > BoostDouble;
typedef BoostDouble InterpolantDouble;

typedef long double long_double;

#ifdef IPM_DOUBLE
    typedef IPM_DOUBLE Double;
#else
    typedef double Double;
#endif

//Change typedef here to use different double type in interior point method.
#ifdef IPM_USE_DOUBLE
typedef Double IPMDouble;
#else
typedef BoostDouble IPMDouble;
#endif

typedef Eigen::Matrix<IPMDouble, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<IPMDouble, Eigen::Dynamic, 1> Vector;

typedef Eigen::Matrix<BoostDouble, Eigen::Dynamic, Eigen::Dynamic>  BoostMatrix;
typedef Eigen::Matrix<BoostDouble, Eigen::Dynamic, 1> BoostVector;

typedef Eigen::Matrix<InterpolantDouble, Eigen::Dynamic, Eigen::Dynamic>  InterpolantMatrix;
typedef Eigen::Matrix<InterpolantDouble, Eigen::Dynamic, 1> InterpolantVector;

typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic>  DoubleMatrix;
typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> DoubleVector;


//TODO: find nicer solution for type casting (does any of the built-in casts might work)

inline DoubleMatrix InterpolantMatrixToMatrix(BoostMatrix &M, DoubleMatrix &){
    DoubleMatrix A(M.rows(), M.cols());
        for (int row = 0; row < M.rows(); ++row) {
            for (int col = 0; col < M.cols(); ++col) {
               A(row,col) = M(row,col).convert_to<Double>();
            } 
        }
    return A;
}

inline BoostMatrix InterpolantMatrixToMatrix(BoostMatrix & M, BoostMatrix &){
    return M;
}

inline DoubleVector InterpolantVectortoVector(BoostVector & v, DoubleVector &){
    DoubleVector w(v.rows());
    for (int i = 0; i < v.rows(); ++i) {
       w(i) = v(i).convert_to<Double>();
    }
    return w;
}

inline BoostVector InterpolantVectortoVector(BoostVector & v, BoostVector &) {
    return v;
}

inline Double InterpolantDoubletoIPMDouble(BoostDouble & d, Double &) {
    return d.convert_to<Double>();
}

inline BoostDouble InterpolantDoubletoIPMDouble(BoostDouble & d, BoostDouble &){
    return d;
}

inline Double InterpolantDoubletoIPMDouble(Double & d, Double &){
    return d;
}

template<class T>
inline T InterpolantDoubletoIPMDouble(Double & d, T){
    return d;
}

inline Vector MatrixToVector(Matrix M){
    assert(M.rows() == M.cols());
    Eigen::Map<Matrix> x(M.data(), M.rows() * M.cols(), 1);
    return x;
}

inline Matrix VectorToSquareMatrix(Vector v, unsigned matrix_dimension){
    assert(v.rows() == matrix_dimension * matrix_dimension);
    Eigen::Map<Matrix> M(v.data(), matrix_dimension, matrix_dimension);
    return M;
}

class Constraints;
class Solution {
public:
    Vector x;
    Vector s;
    IPMDouble centrality;
    IPMDouble gap;
};

//TODO: Need full row rank matrices for IPM. Also, is preprocessing A, e.g. row-echelon form useful?
class Constraints {
public:
    Matrix A;
    Vector b;
    Vector c;

    Constraints (){};
    Constraints(Matrix A_, Vector b_, Vector c_) : A(A_), b(b_), c(c_) {};

    void print(){
        std::cout << "Constraints are as follows. Constraint matrix is A: " << std::endl;
        std::cout << A << std::endl;
        std::cout << "Objective c is " << std::endl;
        std::cout << c.transpose() << std::endl;
        std::cout << "RHS b is " << std::endl;
        std::cout << b.transpose() << std::endl;
    }

    Constraints dual_system(){
        Constraints dual_constraints;
        dual_constraints.c = A.colPivHouseholderQr().solve(b);

        Matrix QR = A.transpose().householderQr().householderQ();
//        std::cout << "Q matrix \n" << QR;

        dual_constraints.A = QR.block(0,A.rows(),QR.rows(), QR.cols() - A.rows()).transpose();
        dual_constraints.b = dual_constraints.A * c;

        //TODO: use different measure to calculate centrality error
        assert((dual_constraints.A * A.transpose()).norm() < 10e-5);
        return dual_constraints;
    }
};

#endif //SOS_UTILS_H




