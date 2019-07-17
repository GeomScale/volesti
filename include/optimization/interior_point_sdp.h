//
// Created by panagiotis on 29/5/2019.
//

#ifndef VOLESTI_INTERIOR_POINT_SDP_H
#define VOLESTI_INTERIOR_POINT_SDP_H


#include "Eigen"
#include "spectrahedron.h"
#include <cmath>
#include <limits>
#include <unsupported/Eigen/MatrixFunctions>


typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<double , Eigen::Dynamic, Eigen::Dynamic> MT;

double inf = std::numeric_limits<double>::infinity();

double f(LMI& lmi, VT& x, double theta) {
    MT A = lmi.evaluate(x);
    double det = A.determinant();

    if (det <= 0)
        return inf;

    double ret = x(x.rows()-1);
    ret -= theta* std::log(det);//A.log().trace();
    return ret;
}

double f(VT& x, double theta, MT& A) {
    double det = A.determinant();

    if (det <= 0)
        return inf;

    double ret = x(x.rows()-1);
    ret -= theta* std::log(det);//A.log().trace();
    return ret;
}


double backtrackingLineSearch(LMI& lmi, const VT &dx, VT& gradient, VT &x, double m) {
    double step_length= 1;
    double alpha = 0.025;
    double beta = 0.7;

    VT eval_x = x + step_length * dx;
//    std::cout << f(lmi, x, m) << " " << f(lmi, eval_x, m)<< "\n";
    double eval_f = f(lmi, x, m);
    double q = alpha * (dx.dot(gradient));

    while (f(lmi, eval_x, m) > eval_f + step_length * q) {
        step_length *= beta;
        eval_x = x + step_length * dx;
    }

    return step_length;
}


void newton(LMI& lmi, VT& x, VT& c, double theta) {
    int dim = x.rows();
    VT gradient(dim);
    MT hessian(dim, dim);
    int steps = 0;

    do {
        MT A = lmi.evaluate(x);
        MT Ainverse = A.inverse();

//        std::cout << lmi.isPositiveDefinite(x)<< " " << lmi.isPositiveSemidefinite(x) << " efe\n";
//        std::cout <<"\n" << x << "\n";

        // compute the gradient
        for (int i = 0; i < dim; i++) {
            const MT& Ai = lmi.getAi(i);
            MT B = Ainverse * Ai;
            gradient(i) = c(i) - theta*B.trace();
        }


        //compute hessian
        for (int i = 0; i < dim ; i++) {
            const MT &Ai = lmi.getAi(i);
            MT AA = Ainverse * Ai * Ainverse;

            for (int j = 0; j < dim; j++) {
                const MT &Aj = lmi.getAi(j);
                MT B = AA * Aj;//TODO save computations
                hessian(i, j) = theta*B.trace();
            }
        }

        VT a = hessian.inverse()*gradient;
        double l = gradient.dot(a);
        if (abs(l) < 0.0000001)
            break;


        VT dx = -1*hessian.inverse() * gradient;

        double step_length = backtrackingLineSearch(lmi, dx, gradient, x, theta);
//        std::cout << step_length << "\n";

        // update point
//        x = x - (1 / (1 + l)) * hessian.inverse() * gradient;
        x = x + step_length * dx;

        // we can also exit if gamma > 0. It suffices for our problem
        if (x(dim-1) < 0) {
//            std::cout << x(dim-1) << "\n";
            break; // it suffices to get gamma above 0
        }
//        std::cout << x(dim-1) << "\n";
    } while (steps++ < 100);
}

void gradientDescent(LMI& lmi, VT& c, VT &x, double m) {
    int dim = x.rows();
    VT gradient(dim);

    do {
        // compute the gradient
        MT A = lmi.evaluate(x);
        MT Ainverse = A.inverse();

//        std::cout << lmi.isPositiveDefinite(x)<< " " << lmi.isPositiveSemidefinite(x) << " efe\n";

//        std::cout << Ainverse << "\n" << x << "\n";

        // compute the gradient
        for (int i = 0; i < dim; i++) {
            const MT& Ai = lmi.getAi(i);
            MT B = Ainverse * Ai;
            gradient(i) = c(i) - m*B.trace();
        }

//        std::cout << A.determinant() << "   " << lmi.isPositiveDefinite(x) << "\n"<< gradient << "-----\n" << x <<"\n";
        double step_length = backtrackingLineSearch(lmi, -1*gradient, gradient, x, m);

//        std::cout << f(x, m, A) << " " << x(dim-1) << " " << step_length << " " << sqrt(gradient.dot(gradient)) <<"\n";

        // we can also exit if s < 0. It suffices for our problem
        if (sqrt(gradient.dot(gradient)) < 0.1) break;

        // compute step length


        // update point
        x -= step_length * gradient;
//        std::cout << f(lmi, x, m) << " " << x(dim-1) << "   " << lmi.isPositiveSemidefinite(x) <<"\n";

        if (x(dim-1) < 0)
            break;
    } while (1);
}

/**
 * Returns a feasible point in the polytope. It uses barrier method to solve the phase I problem to get a feasible point
 *
 * minimize s
 * s.t. f_i(x) <= s
 *
 * The resulting problem is:
 *
 * s - m(log(s - f_1) + log(s - f_2) + ... ), m going to 0

 * If failed at finding a point, throw exception
 *
 * @tparam Point
 * @tparam NT
 * @param polytope
 * @return an feasible point
 */
VT getInteriorPoint(Spectrahedron& spectrahedron) {

    LMI lmi;
    lmi = spectrahedron.getLMI();
    long dim = lmi.getDim();

    // we keep s separate from the other vars, it's easier
    Eigen::EigenSolver<MT> solver;

    MT A0 = lmi.getA0();
    int m = A0.rows();

    solver.compute(A0);
    Eigen::GeneralizedEigenSolver<MT>::ComplexVectorType eivals = solver.eigenvalues();

    double max_gamma = eivals(0).real();

    for (int i = 0; i < eivals.rows(); i++)
        if (eivals(i).real() > max_gamma)
            max_gamma = eivals(i).real();

//    std::cout << max_gamma << "\n";

    if (max_gamma < 0) {
        VT zero;
        zero.setZero(dim);
        return zero;
    }

    MT matrix;
    matrix.setZero(m,m);
    for (int i=0 ; i<m ; i++)
        matrix(i,i) = -1;

    // the initial value for x will be {0, .., 0, max_gamma}
    VT x;
    x.setZero(dim + 1);
    x(dim) = max_gamma + 0.001;
    VT objectiveFunction;
    objectiveFunction.setZero(dim+1);
    objectiveFunction(dim) = 1;

    lmi.addMatrix(matrix);
//    std::cout << lmi.isNegativeDefinite(x)<< " " << lmi.isNegativeSemidefinite(x) << " " <<objectiveFunction<<"efe\n";
    lmi.changeInequalityDirection();
//    std::cout << lmi.isPositiveDefinite(x)<< " " << lmi.isPositiveSemidefinite(x) << " " <<objectiveFunction<<"efe\n";
//    lmi.print();

    int steps = 10;
    double theta = 1;

    while (steps > 0) {
//        std::cout << "========= theta" << theta << "\n" << x << " " << lmi.isPositiveSemidefinite(x)<<"\n";
        if (x(dim) < 0) {
//            std::cout << "========= ending" << theta << "\n" << x << " " << lmi.isPositiveSemidefinite(x)<<"\n";
            break; // it suffices to get gamma above 0
        }
        newton(lmi, x, objectiveFunction, theta);
//        gradientDescent(lmi, objectiveFunction, x, theta);
        steps--;
        theta /= 10;
//        std::cout << lmi.isPositiveSemidefinite(x)<< "efe\n";
    }

    if (x(dim) > 0)
        throw std::string("No internal point was found in the polytope");

    VT _x(dim);
    for (int i=0 ; i<dim ; i++)
        _x(i) = x(i);

    return _x;
}

#endif //VOLESTI_INTERIOR_POINT_SDP_H
