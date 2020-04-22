//
// Created by panagiotis on 29/5/2019.
//

#ifndef VOLESTI_INTERIOR_POINT_H
#define VOLESTI_INTERIOR_POINT_H


#include "Eigen"
#include "polytopes.h"
#include <cmath>
#include <limits>


template<typename NT>
double f(NT s, Eigen::Matrix<NT, Eigen::Dynamic, 1>& x, Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>& A, Eigen::Matrix<NT, Eigen::Dynamic, 1>& b, double m) {
    NT ret = s;
    NT temp;

    for (int i=0 ; i<A.rows(); i++) {
        temp = s - A.row(i).dot(x) + b(i);
        if (temp <= 0)
            return std::numeric_limits<double>::infinity();
        ret -=  m*log(temp);
    }
    return ret;
}

template<typename NT>
double backtrackingLineSearch(NT s, const VT &grad, NT s_grad, MT &A, VT &b, VT &x, double m) {
    double step_length= 1;
    double alpha = 0.025;
    double beta = 0.7;

    Eigen::Matrix<NT, -1, 1> eval = x - step_length * grad;
    NT eval_s = s - step_length*s_grad;

    while (f(eval_s, eval, A, b, m) > f(s, x, A, b, m) - alpha * step_length * (s_grad*s_grad + grad.dot(grad))) {

        step_length *= beta;
        eval = x - step_length * grad;
        eval_s = s - step_length*s_grad;
    }

    return step_length;
}

template<class Point, typename NT>
void gradientDescent(HPolytope<Point> &polytope, MT &A, VT &b, NT &s, VT &x, double m) {
    Eigen::Matrix<NT, -1, 1> denominators;
    denominators.resize(A.rows()); // the denominators from the derivatives of log are common at each iteration

    Eigen::Matrix<NT, -1, 1> grad;
    NT s_grad;
    long dim = polytope.get_mat().cols();

    do {
        // compute the gradient
        for (int i = 0; i < denominators.rows(); i++) {
            denominators(i) =  m / (s - A.row(i).dot(x) + b(i));
        }

        s_grad = 1;
        for (int i = 0; i < denominators.rows(); i++)
            s_grad -=  denominators(i);
        grad.setZero(dim);


        for (int j = 0; j < dim; j++) {
            for (int i = 0; i < A.rows(); i++)
                grad(j) += A(i, j) * denominators(i);
        }

        // we can also exit if s < 0. It suffices for our problem
        if (sqrt(s_grad * s_grad + grad.dot(grad)) < 0.001 || s < 0) break;

        // compute step length
        NT step_length = backtrackingLineSearch(s, grad, s_grad, A, b, x, m);


        // update point
        x -= step_length * grad;
        s -= step_length * s_grad;

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
template<class Point, typename NT>
Point getInteriorPoint(HPolytope<Point>& polytope) {
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;


    MT A = polytope.get_mat();
    VT b = polytope.get_vec();
    long dim = polytope.get_mat().cols();

    // we keep s separate from the other vars, it's easier
    NT s;
    VT x;
    x.resize(dim);

    // the initial value for x will be 1 and for s = max{f_i(x)} + 1
    for (int j=0 ; j<dim ; j++)
        x(j) = 1;

    s = A.row(0).dot(x) - b(0);

    for (int i=1 ; i<A.rows() ; i++) {
        NT _s = A.row(i).dot(x) - b(i);
        if (_s > s)
            s = _s;
    }
    s += 1;

    int steps = 10;
    double m = 1;

    while (steps > 0) {
        gradientDescent(polytope, A, b, s, x, m);
        if (s <0) break; // it suffices to get s below 0
        steps--;
        m /= 10;
    }

    if (s>0)
        throw std::string("No internal point was found in the polytope");

    return Point(x);
}

#endif //VOLESTI_INTERIOR_POINT_H
