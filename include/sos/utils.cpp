#include <Eigen/Dense>
#include <iostream>

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

typedef double Double;

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
    Double centrality;
    Double gap;
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
        std::cout << "Q matrix \n" << QR;

        dual_constraints.A = QR.block(0,A.rows(),QR.rows(), QR.cols() - A.rows()).transpose();
        dual_constraints.b = dual_constraints.A * c;

        //TODO: use different measure to calculate centrality error
        assert((dual_constraints.A * A.transpose()).norm() < 10e-5);
        return dual_constraints;
    }
};



