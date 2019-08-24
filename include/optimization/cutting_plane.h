//
// Created by panagiotis on 24/5/2019.
//

#ifndef VOLESTI_CUTTING_PLANE_H
#define VOLESTI_CUTTING_PLANE_H


#include "polytopes.h"
#include "Eigen"
#include <list>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
#include "samplers.h"
#include "interior_point.h"
#include <queue>
#include "heuristics.h"

int STEPS;

namespace optimization {


    const double ZERO = 0.000000000001;

    /**
     * For the Eigen library
     */
    typedef double NT_MATRIX;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, 1> VT;


    /**
     * Computes the relative error
     * @param approx
     * @param exact
     * @return
     */
    double relative_error(double approx, double exact) {
        return abs((exact - approx) / exact);
    }


    /**
     * A class to save successive approximation values
     */
    class SlidingWindow {
    public:
        std::queue<double> evals;
        int evalsSize; // how many values to store
        int count;
        double min;

        SlidingWindow(int evalSize) {this->evalsSize = evalSize;count = 0;}

        void push(double eval) {
            if (count >= evalsSize) {
                evals.pop();
            }
            else
                count++;

            evals.push(eval);
            min = eval;
        }

        double getRelativeError() {
            if (count < evalsSize)
                return 1;

            return abs((min - evals.front()) / min);
        }

        void half() {
            for (int i=0 ; i<evals.size()/2 ; i++) {
                evals.pop();
                count--;
            }
        }
    };

    /**
     * Returns the normal of a vector
     *
     * @param vector
     * @return
     */
    double getNormal(VT& vector) {
        double normal = 0;
        for (unsigned int i=0; i<vector.rows(); i++)
            normal += vector(i) * vector(i);

        normal=std::sqrt(normal);

        return normal;
    }


    /**
     * Normalizes the given vector
     *
     * @param vector
     * @return the norm of the input vector
     */
    double normalizeVector(VT& vector) {
        double normal = getNormal(vector);

        for (unsigned int i=0; i<vector.rows(); i++)
            vector(i) = vector(i) / normal;

        return normal;
    }


    /**
     * Compute the cutting plane c (x - point) <= 0  ==>  cx <= c point
     * and add it as a new normalized constraint in the existing polytope, in place of its last one, which will be redundant.
     *
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param c a vector holding the coefficients of the object function
     * @param polytope An instance of Hpolytope
     * @param point Where to create the cutting plane
     */
    template<class Point, typename NT>
    void cutPolytope(VT &c, HPolytope<Point> &polytope, Point &point) {

        unsigned int dim = polytope.dimension();
        VT normalizedC = c;
        double normal = normalizeVector(normalizedC);

        //add cx in last row of A
        long j, i;

        for (j = 0, i = polytope.num_of_hyperplanes() - 1 ; j < dim; j++) {
            polytope.put_mat_coeff(i, j, normalizedC(j));
        }

        //add  < c,  point >  in last row of b
        NT _b = c.dot(point.getCoefficients());
        polytope.put_vec_coeff(polytope.num_of_hyperplanes() - 1, _b/normal);
    }


    /**
     * Store in dotProducts the dot product of c with all points in points
     *
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param c a vector holding the coefficients of the object function
     * @param points a list of points
     * @param dotProducts a list to store the dot products
     */
    template<class Point, typename NT>
    void dotProducts(VT &c, std::list<Point> &points, std::list<NT> &dotProducts) {
        for (auto p : points)
            dotProducts.push_back(c.dot(p.getCoefficients()));
    }


    /**
     * Find the point that minimizes the object function out of a list of points
     *
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param c a vector holding the coefficients of the object function
     * @param randPoints a list of points
     * @return point wich minimizes the object function
     */
    template<class Point, typename NT>
    Point getMinimizingPoint(VT &c, std::list<Point> &randPoints) {
        typename std::list<Point>::iterator it = randPoints.begin();

        NT temp, min;
        Point minPoint = *it;

        min = c.dot(it->getCoefficients());
        it++;

        for (; it != randPoints.end(); it++) {
            temp = c.dot(it->getCoefficients());

            if (temp < min) {
                min = temp;
                minPoint = *it;
            }
        }

        return minPoint;
    }


    /**
     * Find the two point that minimize the object function the most out of a list of points
     *
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param c a vector holding the coefficients of the object function
     * @param randPoints a list of points
     * @return a pair of points which minimize the object function (first point minimizes the most)
     */
    template<class Point, typename NT>
    std::pair<Point, Point> getPairMinimizingPoint(VT &c, std::list<Point> &randPoints) {
        class std::list<Point>::iterator points_it = randPoints.begin();
        class std::list<Point>::iterator minPoint, minPoint2;

        NT temp, min, min2;
        minPoint = points_it;

        min = c.dot(points_it->getCoefficients());

        points_it++;

        min2 = c.dot(points_it->getCoefficients());
        minPoint2 = points_it;

        if (min2 < min) {
            temp = min2;
            min2 = min;
            min = temp;
            typename std::list<Point>::iterator tempPoint = minPoint;
            minPoint = minPoint2;
            minPoint2 = tempPoint;
        }

        points_it++;

        for (; points_it != randPoints.end(); points_it++) {
            temp = c.dot(points_it->getCoefficients());

            if (temp < min2) {
                if (temp > min) {
                    min2 = temp;
                    minPoint2 = points_it;
                } else {
                    min2 = min;
                    min = temp;
                    minPoint2 = minPoint;
                    minPoint = points_it;
                }
            }
        }

        return std::pair<Point, Point>(*minPoint, *minPoint2);
    }


    /**
     * Adds one more hyperplane in the polytope. The values of the new polytope are not initialized.
     *
     * @tparam Point class Point
     * @tparam NT the numeric type
     * @param polytope an instance of class HPolytope
     */
    template<class Point, typename NT>
    void addRowInPolytope(HPolytope<Point> &polytope) {
        typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
        typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

        MT A = polytope.get_mat();
        VT b = polytope.get_vec();
        unsigned int dim = polytope.dimension();

        A.conservativeResize(A.rows() + 1, Eigen::NoChange);
        b.conservativeResize(b.rows() + 1);

        polytope.init(dim, A, b);
    }


    /**
     * Prints a message, if verbose is true
     *
     * @param verbose whether or not to print
     * @param msg a message to print
     */
    void print(bool verbose, const char *msg) {
        if (verbose)
            std::cout << msg << std::endl;
    }


    /**
     * Returns the arithmetic mean of points only in the new region of the polytope after cutting it at point p
     *
     * @tparam Point Class Point
     * @tparam NT The numeric type
     * @param c c a vector holding the coefficients of the object function
     * @param points A collection of points
     * @param p the point at which we cut the polytope
     */
    template<class Point, typename NT>
    Point getArithmeticMean(VT &c, std::list<Point> &points, Point& p) {
        class std::list<Point>::iterator points_it = points.begin();

        NT b = c.dot(p.getCoefficients());

        Point point(points_it->dimension());

        int i = 0;

        for (; points_it != points.end(); points_it++)
            if (c.dot(points_it->getCoefficients()) <= b) {
                point = point + *points_it;
                i++;
            }

        point = point  / (double) i;

        return point;
    }


    /**
     * Returns the arithmetic mean of points only in the new region of the polytope after cutting it at point p
     *
     * @tparam Point Class Point
     * @tparam NT The numeric type
     * @param c c a vector holding the coefficients of the object function
     * @param points A collection of points
     * @param p the point at which we cut the polytope
     * @param dotProducts a list with the dot products of c with all elements of points
     */
    template<class Point, typename NT>
    Point getArithmeticMean(VT &c, std::list<Point> &points, Point& p, std::list<NT> &dotProducts) {
        class std::list<Point>::iterator points_it = points.begin();
        typename std::list<NT>::iterator products_it = dotProducts.begin();

        NT b = c.dot(p.getCoefficients());

        Point point(points_it->dimension());

        int i = 0;

        for (; points_it != points.end(); points_it++, products_it++)
            if (*products_it <= b) {
                point = point + *points_it;
                i++;
            }

        point = point / (double) i;

        return point;
    }


    /**
     * Returns the arithmetic mean of points only in the new region of the polytope after cutting it at point p
     *
     * @tparam Point Class Point
     * @tparam NT The numeric type
     * @param c c a vector holding the coefficients of the object function
     * @param points A collection of points
     * @param p the point at which we cut the polytope
     * @param dotProducts a list with the dot products of c with all elements of points
     * @param sum return the sum of all vectors in points
     */
    template<class Point, typename NT>
    Point getArithmeticMean(VT &c, std::list<Point> &points, Point& p, std::list<NT> &dotProducts, VT &sum) {
        class std::list<Point>::iterator points_it = points.begin();
        typename std::list<NT>::iterator products_it = dotProducts.begin();

        NT b = c.dot(p.getCoefficients());

        VT sum_PointsInNewRegion;
        sum_PointsInNewRegion.setZero(points_it->dimension());
        VT sum_PointsInCutRegion;
        sum_PointsInCutRegion.setZero(points_it->dimension());

        int i = 0;

        for (; points_it != points.end(); points_it++, products_it++) {
            if (*products_it <= b) {
                sum_PointsInNewRegion = sum_PointsInNewRegion + points_it->getCoefficients();
                i++;
            }
            else {
                sum_PointsInCutRegion = sum_PointsInCutRegion + points_it->getCoefficients();
            }
        }


        sum = sum_PointsInCutRegion + sum_PointsInNewRegion;
        sum_PointsInNewRegion /= (double) i;
        
        Point point(sum_PointsInNewRegion);

        return point;
    }


    /**
    * Computes the softmax function (gradient of LogSumExp)
    *
    * @param input the input vector
    * @param gradient the returned vector
    * @param weight e = e^weight
    */
    void softMax(VT& input, VT& gradient, double weight) {
        int dim = input.rows();
        std::vector<double> exps(dim, 0);
        double denominator = 0;

        double max = input(0);
        for (int i=1 ; i<dim ; i++)
            if (input(i) > max)
                max = input(i);

        for (int i=0 ; i<dim ; i++) {
            exps[i] = std::exp((input(i) - max)* weight);
            denominator += exps[i];
        }

        gradient.setZero(dim);

        for (int i=0 ; i<dim ; i++) {
            gradient(i) = exps[i] / denominator;
        }
    }



    /**
     * The LogSumExp function
     *
     * @param input
     * @return
     */
    double LogSumExp(VT& input) {
        int dim = input.rows();
        double sum = 0;

        double max = input(0);
        for (int i=1 ; i<dim ; i++)
            if (input(i) > max)
                max = input(i);

        for (int i=0 ; i<dim ; i++) {
            sum += std::exp(input(i) - max);
        }

        return max + log(sum);
    }


    /**
     * Returns the index of the maximum element in the given vector
     * @param vec
     * @return
     */
    int maxElementIndex(VT& vec) {
        int index = 0;
        double max = vec(0);

        for (int i=1 ; i<vec.rows() ; i++)
            if (vec(i) > max) {
                index = i;
                max = vec(i);
            }

        return index;
    }


    /**
     * Returns the index of the minimum element in the given vector
     * @param vec
     * @return
     */
    int minElementIndex(VT& vec) {
        int index = 0;
        double min = vec(0);

        for (int i=1 ; i<vec.rows() ; i++)
            if (vec(i) < min) {
                index = i;
                min = vec(i);
            }

        return index;
    }



    /**
     * Solves the 2 dimensional LP problem
     *
     * max 0x + y
     * subject to a*x + y <= b
     *            x_min <= x <= x_max
     *
     * It's based on Seidel's algorithm for LP problems. The problems I solve always have a feasible solution and y is positive,
     * so the implementation takes these into account to make it more efficient.
     *
     * @param a
     * @param b
     * @param x_min
     * @param x_max
     * @return the coordinates (x,y) of the optimal solution
     */
 //bug
//    std::pair<double, double> Seidel2dimLP(VT &a,  VT &b, double x_min, double x_max, int indexMin, int indexMax) {
//
//
//        // get the intersection points of a constraint with x = x_min, x_max
//        typedef std::pair<double, double> point2;
//
//        point2 max, p, temp, leftmost, rightmost;
//        int chosen = indexMin;
//
//        if (indexMin < indexMax)
//            max = Seidel2dimLP(a, b, x_min, x_max, chosen+1, indexMax);
//        else {
//            rightmost.first = x_min;
//            rightmost.second = - a(indexMax) * x_min + b(indexMax);
//
//            leftmost.first = x_max;
//            leftmost.second = - a(indexMax) * x_max + b(indexMax);
//
//            return rightmost.second < leftmost.second ? leftmost : rightmost;
//        }
//
//        if (a(chosen) * max.first + max.second - b(chosen) <= ZERO)
//            return max;
//
////        p.second = - a(0)*x_min + b(0);
//
//        // coordinates of highest point
////        temp.second = - a(0)*x_max + b(0);
////        max.first = temp.second > p.second ? x_max : x_min;
//
//        // pick the highest y
////        max.second = temp.second > p.second ? temp.second : p.second;
//
//
////        for (int i = 0; i<chosen ; i++) {
////
////             check if constraint is satisfied
////            if (a(i) * max.first + max.second - b(i) <= ZERO)
////                continue;
//
//        // must compute left-most and right-most intersection points with other constraints
//        rightmost.first = x_min;
////            rightmost.second = - a(chosen) * x_min + b(chosen);
//
//        leftmost.first = x_max;
////            leftmost.second = - a(chosen) * x_max + b(chosen);
//
//        for (int j=chosen+1 ; j<indexMax ; j++) {
////                if (i == j) continue;
//
//            if (a(j) - a(chosen) > 0)
//                leftmost.first = (b(j) - b(chosen)) / (a(j) - a(chosen));
//            else
//                rightmost.first = (b(j) - b(chosen)) / (a(j) - a(chosen));
//
////                 must be within boundaries
////                if (p.first < x_min || p.first > x_max) continue;
//
////                p.second = b(chosen) - a(chosen) * p.first;
//
//            // check if constraint j faces right at constraint i
////                temp.first = p.first + 10;
////                temp.second = b(chosen) - a(chosen) * temp.first;
//
////                if (a(j) * temp.first + temp.second - b(j) <= ZERO) {
////                     keep it if it is the most right
////                    if (p.first > rightmost.first) {
////                        rightmost = p;
////                    }
////                } else {
//            // keep it if it is the most left
////                    if (p.first < leftmost.first) {
////                        leftmost = p;
////                    }
////                }
//
//        }
//
//        if (leftmost.first < rightmost.first) { //TODO redundant i constaint?
////                a(i) = a(m-1);
////                b(i) = b(m-1);
////                a.conservativeResize(m-1);
////                b.conservativeResize(m-1);
////                return Seidel2dimLP(a,  b, x_min, x_max, 0, );
////                break;
//            std::cout << "rr\n";
//        }
//
//        // get new highest point
////            if (max.second < leftmost.second) {
////                leftmost.second = b(chosen) - a(chosen) * leftmost.first;
////                max = leftmost;
////            }
////
////            if (max.second < rightmost.second) {
////                rightmost.second = b(chosen) - a(chosen) * rightmost.first;
////                max = rightmost;
////            }
//        if (rightmost.second < leftmost.second) {
//            leftmost.second = b(chosen) - a(chosen) * leftmost.first;
//            max = leftmost;
//        }
//        else {
//            rightmost.second = b(chosen) - a(chosen) * rightmost.first;
//            max = rightmost;
//        }
//
//
////        }
//
//        return max;
//    }



    /**
     * Solves the 2 dimensional LP problem
     *
     * max 0x + y
     * subject to a*x + y <= b
     *            x_min <= x <= x_max
     *
     * @tparam NT
     * @param A
     * @param b
     * @param x
     * @param direction
     * @param lambda_max
     * @param lambda_min
     * @param negativeDistances
     */
    template <typename NT>
    void minimizeAndUpdate(MT &A, const VT &b,  VT &x, const VT &direction, NT lambda_max, NT lambda_min,
                           VT &negativeDistances) {
        VT _b = -1 * negativeDistances;
        VT _a = A*direction;


        double discretize = 10000;
        double stepSize = (lambda_max - lambda_min) / discretize;
        VT step = A * (stepSize * direction);

        double minimizingStep = 1;
        VT t = A * (x+(lambda_min+stepSize)*direction) - b;// = negativeDistances + (A * (lambda_min * direction));
        double min = t(maxElementIndex(t));
        negativeDistances = t;

        for (double i=2 ; i<discretize ; i++) {
            t = A * (x+((lambda_min +  i *stepSize)*direction)) - b;//+= step;
            double n = t(maxElementIndex(t));
//            std::cout << "D \t" << (lambda_min + i *stepSize) << "\t" << n << "\n";
            if (n < min) {
                min = n;
                negativeDistances = t;
                minimizingStep = i;
            } else {
//                std::cout << "D \t" << (lambda_min + i *stepSize) << "\t" << n << "\n";
//                break;
            }
        }
        ///////////////////////////////////////////////////////////////
//
        //        VT v = A * direction;
//        double left = lambda_min;
//        double right = lambda_max;
//        double middle = (right + left) / 2;
//        double min, previous_min;
//
//        VT t = negativeDistances - middle * v;
//        previous_min = t(maxElementIndex(t));
//
//        do {
//
//        } while (relative_error(min, previous_min) < 0.0001);

//        x = x + (lambda_min + minimizingStep * stepSize) * direction;
//        std::cout << lambda_min << "\t"  << lambda_max  << "\t"<< lambda_min + minimizingStep * stepSize << "\t" <<  min << "\n";
//        std::cout << direction.transpose() << "\n";

//        VT _a(4);
//        _a<< 3, 4, 2, -1;
//        VT _b(4);
//        _b << 60, 72, 100, 5;
//        lambda_min = 0;
//        lambda_max = 2;
//        std::pair<double, double> res = Seidel2dimLP(_a, _b, lambda_min, lambda_max,0,_b.rows()-1);


//        x = x + res.first * direction;
//        VT q = -1 * _a * res.first;
//        q = q + _b;
//        std::cout << "\t" << res.first << "\t" <<  res.second  <<"\t" <<q(minElementIndex(q))<< "\n";
/////////////////////////////////////////////////////////////////////////

        lprec *_lp;
        unsigned int dim = 2;

        REAL row[1 + dim]; /* must be 1 more then number of columns ! */

        /* Create a new LP model */
        _lp = make_lp(0, dim);


//        for (int j=1 ; j<=dim ; j++)
//            row[j] = lp.objectiveFunction(j-1); //j must start at 1
        row[1] = 0;
        row[2] = 1;

        set_obj_fn(_lp, row);
        set_verbose(_lp, 2);
        set_add_rowmode(_lp, TRUE);

        for (int i=0 ; i<_a.rows() ; i++) {
//            for (int j=1 ; j<=dim ; j++)
            row[1] = _a(i); //j must start at 1
            row[2] = 1;

            add_constraint(_lp, row, LE, _b(i)); /* constructs the row: +v_1 +2 v_2 >= 3 */
        }

        row[1] = 1; //j must start at 1
        row[2] = 0;

        add_constraint(_lp, row, LE, lambda_max); /* constructs the row: +v_1 +2 v_2 >= 3 */
        row[1] = 1; //j must start at 1
        row[2] = 0;

        add_constraint(_lp, row, GE, lambda_min); /* constructs the row: +v_1 +2 v_2 >= 3 */

        set_add_rowmode(_lp, FALSE);
        set_maxim(_lp);

        for (int j=1 ; j<=dim ; j++)
            set_bounds(_lp, j, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max());

        solve(_lp);
        NT ret = get_objective(_lp);

        REAL solution[dim];
        get_variables(_lp, solution);
//        std::cout << "Optimal solution: " << std::endl;
//        for (int i=0; i < dim; i++)
//            std::cout << "x[" << i+1 << "] = " << solution[i] << std::endl;

        delete_lp(_lp);
        x = x + solution[0] * direction;
        /////////////////////////
//        double min = negativeDistances(maxElementIndex(negativeDistances));
//        diff += (lambda_min - lambda_max) / 25;

//        while (diff > lambda_min) {

//            for (int i = 0; i <  b.rows(); i++)
//                negativeDistances(i) = A.row(i).dot(x+ diff * direction) - b(i);

//                    p = Point((point.getCoefficients() + diff * direction));
//            std::cout << "D \t" << negativeDistances(maxElementIndex(negativeDistances)) << "\t" << diff << "\n";

//            if (negativeDistances(maxElementIndex(negativeDistances)) < min) {
//                min = negativeDistances(maxElementIndex(negativeDistances));
//                best = diff;
//            }
//            diff += (lambda_min - lambda_max) / 25;
//        }
//
//        diff = best;
    }


    /**
     * Computes the euclidean distance between two points.
     *
     * @param v1
     * @param v2
     * @return
     */
    double euclideanDistance(VT& v1, VT& v2) {
        double sum = 0;

        for (int i=0 ; i<v1.rows() ; i++)
            sum += (v1(i) - v2(i))*(v1(i) - v2(i));

        return sqrt(sum);
    }


    template <class Point>
    double euclideanDistance(Point& v1, Point& v2) {
        double sum = 0;

        for (int i=0 ; i<v1.dimension() ; i++)
            sum += (v1[i] - v2[i])*(v1[i] - v2[i]);

        return sqrt(sum);
    }

    /**
     * Normalize the matrix of the polytope
     *
     * @tparam Polytope
     * @param polytope
     */
    template <class Polytope>
    void normalizePolytope(Polytope& polytope) {
        MT A;
        VT b;

        int m;
        b = polytope.get_vec();
        m = b.rows();
        A = polytope.get_mat();

        // normalize
        for (int i = 0; i < m; i++) {
            VT a = A.row(i);
            double norm = getNormal(a);
            A.row(i) /= norm;
            b(i) /= norm;
        }


        polytope.set_mat(A);
        polytope.set_vec(b);
    }

    /**
     * The idea is to move towards the center of the Chebyshev sphere.
     *
     * While at point p, compute vector d = A*p - b
     * Compute the vector s = A^T * softmax(d)
     * Move towards s, that is next point is p' = p + k * s
     * pick k so that you minimize max{a_i*p' - b_i}
     *
     * @tparam Polytope
     * @tparam Point
     * @param polytope
     * @param point
     * @param steps
     */
    template<class Polytope, class Point>
    void escapeStep_ChebyshevCenter(Polytope &polytope, Point &point, int steps) {
        typedef typename Point::FT NT;

//        MT _A = polytope.get_mat();
//        VT _b = polytope.get_vec();

        MT A;
        VT b;
        VT x = point.getCoefficients();
        VT negativeDistances(polytope.get_vec().rows());
        VT gradient;


        int m;
        b = VT(polytope.get_vec());
        m = b.rows();
        A = MT(polytope.get_mat());

        /*****************      INITIALIZE    *****************/


        SlidingWindow slidingWindow(steps/10);

        /*****************      COMPUTE    *****************/

//        VT center = polytope.ComputeInnerBall().first.getCoefficients();
//        VT dd = A*center - b;
//        double radius = dd(maxElementIndex(dd));

        double previous_min = 1;
        double weight = 100;
        bool once = false;

        for (int step = 1; step <= steps; step++) {


            negativeDistances = A * x - b;

//            std::cout << "DISTANCES " << step << "\t Facet: no" << maxElementIndex(negativeDistances) << "\t dist:"
//                      << negativeDistances(maxElementIndex(negativeDistances)) << "\t weight " << weight << "\t"
//                      << euclideanDistance(center, x) << "\t" <<radius<< "\n";



            // find direction
            softMax(negativeDistances, gradient, weight);
            VT direction = A.transpose() * gradient;

            // how far can we travel
            std::pair<NT, NT> dbpair = polytope.line_intersect(Point(x), Point(direction));
            NT lambda_max = dbpair.first;
            NT lambda_min = dbpair.second;


//            VT _b = -1 * negativeDistances;
//            VT _a = A*direction;
//            std::pair<double, double> res = Seidel2dimLP(_a, _b, lambda_min, lambda_max);
//            x = x + res.first * direction;
            minimizeAndUpdate(A,b,  x, direction, lambda_max, lambda_min, negativeDistances);
            double res = negativeDistances(maxElementIndex(negativeDistances));

            if (previous_min == 1)
                previous_min = res;
            else {
                if (relative_error(previous_min, res) < 0.1) {
                    weight *= 10;

                    if (weight > 10000) {
                        weight = 10;

//                        if (once) {
//                        std::cout << A.row(A.rows() - 1) << "\n";
//                        std::cout << A.row(maxElementIndex(negativeDistances)) << "\n";
//                            VT maxFacet = A.row(maxElementIndex(negativeDistances));
//
//                            double sum = 0;
//                            int nonzero = -1;
//                            int no = 0;
//
//                            for (int i = 0; i < A.cols(); i++) {
//
//                                if (maxFacet(i) != 0 && nonzero == -1) {
//                                    nonzero = i;
//                                    continue;
//                                }
//
//                                if (maxFacet(i) != 0) {
//                                    sum += maxFacet(i) * 2;
//                                    no++;
//                                }
//                            }
//
//                            if (no != 0) {
//                                A.conservativeResize(A.rows() + 1, A.cols());
//
//
//                                A(A.rows() - 1, nonzero) = sum / no;
//                                b.conservativeResize(b.rows() + 1);
//                                b(b.rows() - 1) = A.row(A.rows() - 1).dot(x) + 0.0001;
//
//                                m = b.rows();
//
//                                for (int i = 0; i < m; i++) {
//                                    VT a = A.row(i);
//                                    double norm = getNormal(a);
//                                    A.row(i) /= norm;
//                                    b(i) /= norm;
//                                }
//                            std::cout << A.row(A.rows() - 1) << "\n" << maxFacet.transpose() << "\n"
//                                      << A.row(A.rows() - 1).dot(maxFacet) << "\n" << A.row(maxElementIndex(negativeDistances))
//                                      << "\n";
//                                std::cout <<  "\n\n\n";
//                                polytope.set_mat(A);
//                                polytope.set_vec(b);
//                                once = false;
//                            }
//                        }

                    }


                }

                previous_min = res;

            }
            slidingWindow.push(res);

            if (slidingWindow.getRelativeError() < 0.01)
                break;
        }

        point = Point(x);

//        assert(false);

//        polytope.set_mat(_A);
//        polytope.set_vec(_b);

//        assert(polytope.is_in(Point(x)));

    }



    /**
     * Use billard walk as an escape step. Set the point to be the center of the longest trajectory.
     *
     * @tparam Polytope
     * @tparam Point
     * @param polytope
     * @param point
     * @param walkLength
     */
    template<class Polytope, class Point>
    void escapeStep_BilliardWalk(Polytope& polytope, Point& point, int walkLength) {
        typedef typename Point::FT NT;

        MT A = polytope.get_mat();
        VT b = polytope.get_vec();
        VT x = point.getCoefficients();

        VT negativeDistances = A * x - b;
        VT gradient;

        // find direction
        softMax(negativeDistances, gradient, 10000);
        VT direction = A.transpose() * gradient;

        // the center and the length of each trajectory
        std::vector<std::pair<VT, double> > lengths;

        for (int i=0 ; i<walkLength ; i++) {

            // how far can we travel
            Point facet_first, facet_second;
            std::pair<NT, NT> dbpair = polytope.line_intersect(Point(x), Point(direction), facet_first, facet_second);

            double lambda = abs(dbpair.first) > abs(dbpair.second) ? dbpair.first : dbpair.second;
            VT s = abs(dbpair.first) > abs(dbpair.second) ? facet_first.getCoefficients() : facet_second.getCoefficients();

            VT endPoint1 = x + dbpair.first * direction;
            VT endPoint2 = x + dbpair.second * direction;
            double length = euclideanDistance(endPoint1, endPoint2);
            VT center = (endPoint1 + endPoint2)/2;

//            std::cout << "length = " << length << ", is in = " << polytope.is_in(Point(center)) << " lambdas:" << dbpair.first << ", " << dbpair.second << " " <<
//                   facet_first.getCoefficients()(0) << " " << facet_second.getCoefficients()(0) << "\n";
            lengths.push_back(std::pair<VT, double>(center, length));

            x = x + 0.9*lambda* direction;
            direction = direction - 2* direction.dot(s) * s;
        }

        double max = lengths[0].second;
        int max_index = 0;

        for (int i=0 ; i<lengths.size() ; i++){
            if (lengths[i].second > max) {
                max = lengths[i].second;
                max_index = i;
            }
        }

        point = Point(lengths[max_index].first);

    }

    /**
     * Use billard walk as an escape step. Set the point to be the arithmetic mean of points sampled.
     *
     * @tparam Polytope
     * @tparam Point
     * @param polytope
     * @param point
     * @param walkLength
     */
    template<class Polytope, class Point>
    void escapeStep_BilliardWalkCenter(Polytope& polytope, Point& point, int walkLength) {
        typedef typename Point::FT NT;

        MT A = polytope.get_mat();
        VT b = polytope.get_vec();
        VT x = point.getCoefficients();

        VT negativeDistances = A * x - b;
        VT gradient;

        // find direction
        softMax(negativeDistances, gradient, 10000);
        VT direction = A.transpose() * gradient;

        // the center and the length of each trajectory
        std::vector<std::pair<VT, double> > lengths;

        for (int i=0 ; i<walkLength ; i++) {

            // how far can we travel
            Point facet_first, facet_second;
            std::pair<NT, NT> dbpair = polytope.line_intersect(Point(x), Point(direction), facet_first, facet_second);

            double lambda = abs(dbpair.first) > abs(dbpair.second) ? dbpair.first : dbpair.second;
            VT s = abs(dbpair.first) > abs(dbpair.second) ? facet_first.getCoefficients() : facet_second.getCoefficients();

            VT endPoint1 = x + dbpair.first * direction;
            VT endPoint2 = x + dbpair.second * direction;
            double length = euclideanDistance(endPoint1, endPoint2);
            VT center = (endPoint1 + endPoint2)/2;

//            std::cout << "length = " << length << ", is in = " << polytope.is_in(Point(center)) << " lambdas:" << dbpair.first << ", " << dbpair.second << " " <<
//                      facet_first.getCoefficients()(0) << " " << facet_second.getCoefficients()(0) << "\n";
            lengths.push_back(std::pair<VT, double>(center, length));

            x = x + 0.9*lambda* direction;
            direction = direction - 2* direction.dot(s) * s;
        }

        VT sum;
        sum.setZero(direction.rows());

        for (int i=0 ; i<lengths.size() ; i++)
            sum = sum + lengths[i].first;

        point = Point(sum / (double) lengths.size());

    }

    //vectors normalized
    VT getVectorWithAngle(VT& s, double cos) {
        int dim = s.rows();
        VT v(dim);

        int firstNonZero = -1;
        double sum = 0;

        for (int i=0 ; i<dim; i++) {
            v(i) = 1;

            if (s(i) != 0 && firstNonZero == -1) {
                    firstNonZero = i;
                    continue;
            }

            sum += s(i);
        }

        v(firstNonZero) = (cos - sum) / s(firstNonZero);
        return v;
    }

    /**
     * Use billard walk as an escape step. Set the point to be the one that minimizes the objective function.
     *
     * @tparam Polytope
     * @tparam Point
     * @param polytope
     * @param point
     * @param walkLength
     * @param c
     */
    template<class Polytope, class Point>
    void escapeStep_BilliardWalk(Polytope& polytope, Point& point, int walkLength, VT c) {
        typedef typename Point::FT NT;

        MT A = polytope.get_mat();
        VT b = polytope.get_vec();
        VT x = point.getCoefficients();

        VT negativeDistances = A * x - b;
        VT gradient;

        // find direction
        softMax(negativeDistances, gradient, 10000);
        VT direction = c;//A.transpose() * gradient;

        // the center and the length of each trajectory
        std::vector<std::pair<VT, VT> > lengths;

        for (int i=0 ; i<walkLength ; i++) {

            // how far can we travel
            Point facet_first, facet_second;
            std::pair<NT, NT> dbpair = polytope.line_intersect(Point(x), Point(direction), facet_first, facet_second);

//            double lambda = abs(dbpair.first) > abs(dbpair.second) ? dbpair.first : dbpair.second;
            VT s = abs(dbpair.first) > abs(dbpair.second) ? facet_first.getCoefficients() : facet_second.getCoefficients();

            double lambda = (dbpair.first + dbpair.second) / 2;
//            VT endPoint1 = x + dbpair.first * direction;
//            VT endPoint2 = x + dbpair.second * direction;
//            double length = euclideanDistance(endPoint1, endPoint2);
            VT center = x + lambda * direction;

//            std::cout << "length = " << length << ", is in = " << polytope.is_in(Point(center)) << " lambdas:" << dbpair.first << ", " << dbpair.second << " " <<
//                      c.dot(endPoint1) << " " << c.dot(endPoint2) << "\n";

            lengths.push_back(std::pair<VT, VT>(center, s));
//            lengths.push_back(std::pair<VT, VT>(endPoint1, s));
//            lengths.push_back(std::pair<VT, VT>(endPoint2, s));

            x = x + 0.9*lambda* direction;
            direction = direction - 2* direction.dot(s) * s;
        }

        double min = lengths[0].first.dot(c);
        int min_index = 0;

        for (int i=0 ; i<lengths.size() ; i++){
            if (lengths[i].first.dot(c) < min) {
                min = lengths[i].first.dot(c);
                min_index = i;
            }
        }

        point = Point(lengths[min_index].first);
//        std::cout << lengths[min_index].first.dot(c) << "\n";
    }



    /**
     * Check if new point p is a better approximation than min1, min2 and if yes change min1, min2
     *
     *
     * @tparam Point
     * @tparam NT
     * @param p
     * @param dotProduct1 dot product of min1 with objective function
     * @param dotProduct2 dot product of min2 with objective function
     * @param min1 the point that minimizes the objective function
     * @param min2 the second best point
     * @param newProduct dot product of p with objective function
     * @param changedMin1 set true if changed min1
     * @param changedMin2 set true if changed min2
     */
    template <class Point, typename NT>
    void getNewMinimizingPoints(const Point &p, NT& dotProduct1, NT& dotProduct2, Point &min1, Point &min2, NT newProduct, bool& changedMin1, bool& changedMin2) {

        if (newProduct < dotProduct2) {
            if (newProduct < dotProduct1) {
                dotProduct2 = dotProduct1;
                min2 = min1;
                dotProduct1 = newProduct;
                min1 = p;
                changedMin1 = true;
            } else {
                dotProduct2 = newProduct;
                min2 = p;
                changedMin2 = true;
            }
        }

    }

    /**
     * Check if new point p is a better approximation than min1, min2 and if yes change min1, min2
     *
     *
     * @tparam Point
     * @tparam NT
     * @param p
     * @param dotProduct1 dot product of min1 with objective function
     * @param dotProduct2 dot product of min2 with objective function
     * @param min1 the point that minimizes the objective function
     * @param min2 the second best point
     * @param newProduct dot product of p with objective function
     */
    template <class Point, typename NT>
    void getNewMinimizingPoints(const Point &p, NT& dotProduct1, NT& dotProduct2, Point &min1, Point &min2, NT newProduct) {

        if (newProduct < dotProduct2) {
            if (newProduct < dotProduct1) {
                dotProduct2 = dotProduct1;
                min2 = min1;
                dotProduct1 = newProduct;
                min1 = p;
            } else {
                dotProduct2 = newProduct;
                min2 = p;
            }
        }

    }



    /**
     * Generate random points and return the two that minimize the objective function.
     * Also returns in list endPoints the intersection points of the polytope with each ray of the random walk
     *
     * @tparam Polytope
     * @tparam PointList
     * @tparam Parameters
     * @tparam Point
     * @param P
     * @param c c the objective function
     * @param p p a interior point
     * @param rnum # of points to
     * @param rnum num # of points to generate
     * @param endPointEveryKSteps how often to collect endpoints
     * @param endPoints
     * @param var defines which walk to use
     * @return (p1, p2) p1 minimizes c the most and p2 follows
     */
    template<class Polytope, class PointList, class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_boundary(Polytope &P,
                                                              VT& c,
                                                              Point &p,   // a point to start
                                                              unsigned int rnum,
                                                              unsigned int endPointEveryKSteps,
                                                              PointList &endPoints,
                                                              Parameters &var)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;
//        assert(P.is_in(p));

        int dim = p.dimension();
        VT b = P.get_vec();

        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);
        std::vector<NT> lamdas(P.num_of_hyperplanes(), NT(0));
        unsigned int rand_coord, rand_coord_prev;
        NT kapa, ball_rad = var.delta;
        Point p_prev = p;

        if (var.cdhr_walk) {//Compute the first point for the CDHR
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
            p_prev = p;
            p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        } else
            hit_and_run(p, P, var);

        min1 = p;
        NT dotProduct1 = min1.getCoefficients().dot(c);

        if (var.cdhr_walk) {//Compute the first point for the CDHR
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
            p_prev = p;
            p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        } else
            hit_and_run(p, P, var);

        min2 = p;
        NT dotProduct2 = min2.getCoefficients().dot(c);
        NT currentProduct = dotProduct2;

        if (dotProduct1 > dotProduct2) {
            NT temp = dotProduct1;
            dotProduct1 = dotProduct2;
            dotProduct2 = temp;
            Point t = min1;
            min1 = min2;
            min2 = t;
        }

        std::pair<NT, NT> bpair;

        NT diff;
        p1 = p2 = p;
        currentProduct = p.getCoefficients().dot(c);
        int count = 0;

        for (unsigned int i = 1; i <= rnum; ++i) {

            if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas, bpair);
                diff = bpair.first + kapa * (bpair.second - bpair.first);
                currentProduct += c(rand_coord) * diff;
                p1 = p2 = p;
                p1.set_coord(rand_coord, p1[rand_coord] + bpair.first);
                p2.set_coord(rand_coord, p2[rand_coord] + bpair.second);
            } else {
                hit_and_run(p, P, var, p1, p2);
                currentProduct = p.getCoefficients().dot(c);
            }

            count++;
            if (count == endPointEveryKSteps) {
                endPoints.push_back(p1);
                endPoints.push_back(p2);
                count = 0;
            }


            getNewMinimizingPoints(p, dotProduct1, dotProduct2, min1, min2, currentProduct);
//        assert(P.is_in(min1));
//        assert(P.is_in(min2));
//        std::cout << min1.getCoefficients().dot(c) <<"\n";
        }

        p = (min1 + min2)/2;

        return std::pair<Point, Point>(min1, min2);
    }


    /**
     * Generate random points and return the two that minimize the objective function
     *
     * @tparam Polytope
     * @tparam Parameters
     * @tparam Point
     * @param P
     * @param c the objective function
     * @param p a interior point
     * @param rnum # of points to generate
     * @param var defines which walk to use
     * @return (p1, p2) p1 minimizes c the most and p2 follows
     */
    template<class Polytope, class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator(Polytope &P,
                                                     VT& c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);
        std::vector<NT> lamdas(P.num_of_hyperplanes(), NT(0));
        unsigned int rand_coord, rand_coord_prev;
        NT kapa, ball_rad = var.delta;
        Point p_prev = p;

        if (var.cdhr_walk) {//Compute the first point for the CDHR
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
            p_prev = p;
            p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        } else
            hit_and_run(p, P, var);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        if (var.cdhr_walk) {//Compute the first point for the CDHR
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
            p_prev = p;
            p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        } else
            hit_and_run(p, P, var);

        min2 = p;
        NT minProduct2 = min2.getCoefficients().dot(c);
        NT newProduct = minProduct2;

        if (minProduct1 > minProduct2) {
            NT temp = minProduct1;
            minProduct1 = minProduct2;
            minProduct2 = temp;
            Point t = min1;
            min1 = min2;
            min2 = t;
        }

        std::pair<NT, NT> bpair;

        // this point will be the end point of the segment of the minimizing point, that lies in the polytope after the cut
        // we will use it to get an interior point to start the random walk at the next phase
        Point boundaryMin1 = min1;
        Point boundaryMin2 = min2;

        // begin sampling

        for (unsigned int i = 1; i <= rnum ; ++i) {

            // get next point
            if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);

                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas, bpair);
                newProduct += c(rand_coord) * (bpair.first + kapa * (bpair.second - bpair.first));
            } else {
                hit_and_run(p, P, var, p1, p2);
                newProduct = p.getCoefficients().dot(c);
            }


            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

//            if (changedMin1){
//                // if the new point is the new min update the boundary point
//
//                boundaryMin2 = boundaryMin1;
//                if (var.cdhr_walk) {
//                    if (c(rand_coord) > 0) {
//                        boundaryMin1 = p_prev;
//                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.second);
//                    } else {
//                        boundaryMin1 = p_prev;
//                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.first);
//                    }
//
//                } else {
//                    boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//                }
//            }
//            else if (changedMin2) {
//                if (var.cdhr_walk) {
//                    if (c(rand_coord) > 0) {
//                        boundaryMin2 = p_prev;
//                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.second);
//                    } else {
//                        boundaryMin2 = p_prev;
//                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.first);
//                    }
//
//                } else {
//                    boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//                }
//            }

//        assert(P.is_in(min1));
//        assert(P.is_in(min2));
//            std::cout << min1.getCoefficients().dot(c) << "\t" << minProduct1 <<"\n";
        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

//TODO done need boundary perhaps?
        // find an interior point to start the next phase
        Point _p =  min1*0.50;
        Point _p1 = min2*0.50;
//        Point _p2 =  boundaryMin1*0.20;
//        Point _p3 = boundaryMin2*0.20;
//        p = (min1 + min2)/2;
        p = _p + _p1;// + _p2 + _p3;

//        std::cout << min1.getCoefficients().dot(c) << "\t" << P.is_in(p) << "\t" << min2.getCoefficients().dot(c) << "\n";
//        std::cout << min2.getCoefficients().dot(c) <<"\n";
//        std::cout << "b " << boundaryMin1.getCoefficients().dot(c) << " " << P.is_in(boundaryMin1) << "\n";

        return std::pair<Point, Point>(min1, min2);
    }

    /**
     * Generate random points and return the two that minimize the objective function
     *
     * @tparam Polytope
     * @tparam Parameters
     * @tparam Point
     * @param P
     * @param c the objective function
     * @param p a interior point
     * @param rnum # of points to generate
     * @param var defines which walk to use
     * @return (p1, p2) p1 minimizes c the most and p2 follows
     */
    template<class Polytope, class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator(Polytope &P,
                                                     VT& c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                     std::list<Point> points)
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);
        std::vector<NT> lamdas(P.num_of_hyperplanes(), NT(0));
        unsigned int rand_coord, rand_coord_prev;
        NT kapa, ball_rad = var.delta;
        Point p_prev = p;

        if (var.cdhr_walk) {//Compute the first point for the CDHR
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
            p_prev = p;
            p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        } else
            hit_and_run(p, P, var);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        if (var.cdhr_walk) {//Compute the first point for the CDHR
            rand_coord = uidist(rng);
            kapa = urdist(rng);
            std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
            p_prev = p;
            p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));
        } else
            hit_and_run(p, P, var);

        min2 = p;
        NT minProduct2 = min2.getCoefficients().dot(c);
        NT newProduct = minProduct2;

        if (minProduct1 > minProduct2) {
            NT temp = minProduct1;
            minProduct1 = minProduct2;
            minProduct2 = temp;
            Point t = min1;
            min1 = min2;
            min2 = t;
        }

        std::pair<NT, NT> bpair;

        // this point will be the end point of the segment of the minimizing point, that lies in the polytope after the cut
        // we will use it to get an interior point to start the random walk at the next phase
        Point boundaryMin1 = min1;
        Point boundaryMin2 = min2;

        // begin sampling

        for (unsigned int i = 1; i <= rnum ; ++i) {

            // get next point
            if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);

                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas, bpair);
                newProduct += c(rand_coord) * (bpair.first + kapa * (bpair.second - bpair.first));
            } else {
                hit_and_run(p, P, var, p1, p2);
                newProduct = p.getCoefficients().dot(c);
            }
            points.push_back(p);

            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

//            if (changedMin1){
                // if the new point is the new min update the boundary point

//                boundaryMin2 = boundaryMin1;
//                if (var.cdhr_walk) {
//                    if (c(rand_coord) > 0) {
//                        boundaryMin1 = p_prev;
//                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.second);
//                    } else {
//                        boundaryMin1 = p_prev;
//                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.first);
//                    }
//
//                } else {
//                    boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//                }
//            }
//            else if (changedMin2) {
//                if (var.cdhr_walk) {
//                    if (c(rand_coord) > 0) {
//                        boundaryMin2 = p_prev;
//                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.second);
//                    } else {
//                        boundaryMin2 = p_prev;
//                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.first);
//                    }
//
//                } else {
//                    boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//                }
//            }

//        assert(P.is_in(min1));
//        assert(P.is_in(min2));
//            std::cout << min1.getCoefficients().dot(c) << "\t" << minProduct1 <<"\n";
        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

//TODO done need boundary perhaps?
        // find an interior point to start the next phase
        Point _p =  min1*0.50;
        Point _p1 = min2*0.50;
//        Point _p2 =  boundaryMin1*0.20;
//        Point _p3 = boundaryMin2*0.20;
//        p = (min1 + min2)/2;
        p = _p + _p1;// + _p2 + _p3;

//        std::cout << min1.getCoefficients().dot(c) << "\t" << P.is_in(p) << "\t" << min2.getCoefficients().dot(c) << "\n";
//        std::cout << min2.getCoefficients().dot(c) <<"\n";
//        std::cout << "b " << boundaryMin1.getCoefficients().dot(c) << " " << P.is_in(boundaryMin1) << "\n";

        return std::pair<Point, Point>(min1, min2);
    }

    /**
     * Generate random points and return the two that minimize the objective function
     *
     * @tparam Polytope
     * @tparam Parameters
     * @tparam Point
     * @param P
     * @param c the objective function
     * @param p a interior point
     * @param rnum # of points to generate
     * @param var defines which walk to use
     * @return (p1, p2) p1 minimizes c the most and p2 follows
     */
    template<class Polytope, class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator(Polytope &P,
                                                     VT& c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                     std::list<Point> points,
                                                     MT& covarianceMatrix)
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);

        hit_and_run_sampled_covariance_matrix(p, P, var, covarianceMatrix);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        hit_and_run_sampled_covariance_matrix(p, P, var, covarianceMatrix);

        min2 = p;
        NT minProduct2 = min2.getCoefficients().dot(c);
        NT newProduct = minProduct2;

        if (minProduct1 > minProduct2) {
            NT temp = minProduct1;
            minProduct1 = minProduct2;
            minProduct2 = temp;
            Point t = min1;
            min1 = min2;
            min2 = t;
        }

        std::pair<NT, NT> bpair;

        // this point will be the end point of the segment of the minimizing point, that lies in the polytope after the cut
        // we will use it to get an interior point to start the random walk at the next phase
        Point boundaryMin1 = min1;
        Point boundaryMin2 = min2;

        // begin sampling

        for (unsigned int i = 1; i <= rnum ; ++i) {

            // get next point
            hit_and_run_sampled_covariance_matrix(p, P, var, covarianceMatrix);
            newProduct = p.getCoefficients().dot(c);

            points.push_back(p);

            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

//            if (changedMin1){
                // if the new point is the new min update the boundary point
//
//                boundaryMin2 = boundaryMin1;
//                if (var.cdhr_walk) {
//                    if (c(rand_coord) > 0) {
//                        boundaryMin1 = p_prev;
//                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.second);
//                    } else {
//                        boundaryMin1 = p_prev;
//                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.first);
//                    }
//
//                } else {
//                    boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//                }
//            }
//            else if (changedMin2) {
//                if (var.cdhr_walk) {
//                    if (c(rand_coord) > 0) {
//                        boundaryMin2 = p_prev;
//                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.second);
//                    } else {
//                        boundaryMin2 = p_prev;
//                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.first);
//                    }

//                } else {
//                    boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//                }
//            }

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

        // find an interior point to start the next phase
        Point _p =  min1*0.50;
        Point _p1 = min2*0.50;
        p = _p + _p1;

        return std::pair<Point, Point>(min1, min2);
    }


    /**
 * Generate random points and return the two that minimize the objective function
 *
 * @tparam Polytope
 * @tparam Parameters
 * @tparam Point
 * @param P
 * @param c the objective function
 * @param p a interior point
 * @param rnum # of points to generate
 * @param var defines which walk to use
 * @return (p1, p2) p1 minimizes c the most and p2 follows
 */
    template<class Polytope, class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_billiard(Polytope &P,
                                                     VT& c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                      double diameter)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);
        std::vector <NT> lamdas(P.num_of_hyperplanes(), NT(0));
        std::vector <NT> Av(P.num_of_hyperplanes(), NT(0));
        double lambda;

        Point p1(dim), p2(dim), min1(dim), min2(dim);
        unsigned int rand_coord, rand_coord_prev;
        NT kapa, ball_rad = var.delta;
        Point p_prev = p;


        billiard_walk(P, p, diameter, lamdas, Av, lambda, var, true);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);


        billiard_walk(P, p, diameter, lamdas, Av, lambda,  var);

        min2 = p;
        NT minProduct2 = min2.getCoefficients().dot(c);
        NT newProduct = minProduct2;

        if (minProduct1 > minProduct2) {
            NT temp = minProduct1;
            minProduct1 = minProduct2;
            minProduct2 = temp;
            Point t = min1;
            min1 = min2;
            min2 = t;
        }

        std::pair<NT, NT> bpair;

        // begin sampling

        for (unsigned int i = 1; i <= rnum ; ++i) {


            billiard_walk(P, p, diameter, lamdas, Av, lambda,  var);
            newProduct = p.getCoefficients().dot(c);


            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

        // find an interior point to start the next phase
        Point _p =  min1*0.50;
        Point _p1 = min2*0.50;

        p = _p + _p1;// + _p2 + _p3;

        return std::pair<Point, Point>(min1, min2);
    }

    template <class Point, typename NT>
    std::pair<NT, NT>
    getCuttingInteriorPoints(const VT &objectiveFunction, Point &interiorPoint, unsigned int rand_coord,
                             const std::pair<NT, NT> &bpair) {

        NT cutAt, diff;
        if (objectiveFunction(rand_coord) > 0) {

            if (abs(bpair.second) < ZERO)
                return std::pair<NT, NT>(0, 0);

            if (relative_error(interiorPoint[rand_coord] + bpair.second, interiorPoint[rand_coord]) > 0.001) {
                cutAt =  0.7*bpair.second;
                diff = 0.9 * bpair.second;
            }
            else {
                cutAt =  0.3*bpair.second;
                diff = 0.6*bpair.second;
            }
        } else {
            if (abs(bpair.first) < ZERO)
                return std::pair<NT, NT>(0, 0);

            if (relative_error(interiorPoint[rand_coord] + bpair.first, interiorPoint[rand_coord]) > 0.001) {
                cutAt = 0.7*bpair.first;
                diff = 0.9 * bpair.first;
            }
            else {
                cutAt = 0.3*bpair.first;
                diff = 0.6*bpair.first;
            }
        }

        return std::pair<NT, NT>(cutAt, diff);
    }


    /**
     * Solve the linear program
     *
     *      min cx
     *      s.t. Ax <= b
     *
     * In this implementation, at each phase of the algorithm, I generate a new point with the cdhr random walk and cut the polytope
     *
     * @tparam Parameters
     * @tparam Point
     * @tparam NT
     * @param polytope
     * @param objectiveFunction
     * @param parameters
     * @param error
     * @param maxSteps
     * @param initial
     * @return
     */
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method_new(HPolytope<Point> polytope, VT& objectiveFunction, Parameters parameters, const NT error,
                                                  const unsigned int maxSteps, Point& initial) {


        typedef typename Parameters::RNGType RNGType;


        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        unsigned int walk_len = parameters.walk_steps;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 1;
        int dim = initial.dimension();


        SlidingWindow slidingWindow(2*dim);
        std::pair<Point, Point> minimizingPoints;
        bool escape = false;
        bool escape_chebysev = false;
        bool escape_billiard = false;
        Point interiorPoint = initial;


        // initialize the cdhr random walk
        RNGType &rng = parameters.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);
        std::vector<NT> lamdas(polytope.num_of_hyperplanes(), NT(0));
        unsigned int rand_coord, rand_coord_prev;
        NT kapa, ball_rad = parameters.delta;
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        std::pair<NT, NT> bpair = polytope.line_intersect_coord(interiorPoint, rand_coord, lamdas);
        Point interiorPoint_prev = interiorPoint;


        std::pair<NT, NT> res = getCuttingInteriorPoints(objectiveFunction, interiorPoint, rand_coord, bpair);
        double cutAt = res.first;
        NT diff = res.second;

        // cut the polytope
        interiorPoint.set_coord(rand_coord, interiorPoint[rand_coord] + cutAt);

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        cutPolytope<Point, NT>(objectiveFunction, polytope, interiorPoint);
        normalizePolytope(polytope);

        // prepare point for next iteration and update min value
        interiorPoint.set_coord(rand_coord, interiorPoint[rand_coord] - cutAt + diff);

        lamdas = std::vector<double>(polytope.num_of_hyperplanes(), NT(0));
        rand_coord = uidist(rng);
        kapa = urdist(rng);
        bpair = polytope.line_intersect_coord(interiorPoint, rand_coord, lamdas);
        interiorPoint_prev = interiorPoint;
        interiorPoint.set_coord(rand_coord, interiorPoint[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));

        NT dotProduct = interiorPoint.getCoefficients().dot(objectiveFunction);
        slidingWindow.push(dotProduct);


        res = getCuttingInteriorPoints(objectiveFunction, interiorPoint, rand_coord, bpair);
        cutAt = res.first;
        diff = res.second;


        double previous_dotProduct;
        int stepsSinceLastEscape = 0;
        double objectiveFunctionNorm = getNormal(objectiveFunction);
        bool stuck = false;
        int stuck_try_coord;
        bool escapeDidntWork = false;
        bool triedEscaping = false;

        do {
            previous_dotProduct = dotProduct;

            if (escape) {
                if (escape_chebysev)
                    escapeStep_ChebyshevCenter(polytope, interiorPoint, 3);
                else if (escape_billiard)
                    escapeStep_BilliardWalk(polytope, interiorPoint, walk_len, objectiveFunction);

                rand_coord = uidist(rng);
                kapa = urdist(rng);
                bpair = polytope.line_intersect_coord(interiorPoint, rand_coord, lamdas);
                interiorPoint_prev = interiorPoint;
                interiorPoint.set_coord(rand_coord, interiorPoint[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));

                dotProduct = interiorPoint.getCoefficients().dot(objectiveFunction);

                escape = false;
                triedEscaping = true;
                stepsSinceLastEscape = 0;
            }
            else {

                // get next point
                rand_coord_prev = rand_coord;

                if (!stuck)
                    rand_coord = uidist(rng);
                else {
                    rand_coord = stuck_try_coord;
                    stuck_try_coord++;
                }

                kapa = urdist(rng);

                bpair = polytope.line_intersect_coord(interiorPoint, interiorPoint_prev, rand_coord, rand_coord_prev,
                                                      lamdas);
                interiorPoint_prev = interiorPoint;

                res = getCuttingInteriorPoints(objectiveFunction, interiorPoint, rand_coord, bpair);
                cutAt = res.first;
                diff = res.second;

                if (diff != 0 && cutAt != 0) {

                    polytope.put_vec_coeff(polytope.num_of_hyperplanes() - 1, (dotProduct + objectiveFunction(rand_coord) * cutAt) / objectiveFunctionNorm);

                    // prepare point for next iteration and update min value
                    interiorPoint.set_coord(rand_coord, interiorPoint[rand_coord] + diff);
                    dotProduct += objectiveFunction(rand_coord) * diff;

                    //update lambdas for random walk
                    VT b = polytope.get_vec();
                    MT A = polytope.get_mat();
                    int at = b.rows() - 1;
                    NT sum_nom = b(at);
                    NT sum_denom = A(at, rand_coord);
                    VT r_coeffs = interiorPoint_prev.getCoefficients();
                    for (int j = 0; j < interiorPoint_prev.dimension(); j++) {
                        sum_nom -= A(at, j) * r_coeffs(j);
                    }
                    lamdas[at] = sum_nom;

                }
            }

            slidingWindow.push(dotProduct);

            if (slidingWindow.getRelativeError() < error) {
                if (triedEscaping) {
                    if (!stuck) {
                        stuck = true;
                        stuck_try_coord = 0;
                        escape = false;
                    }
                    else if (stuck_try_coord == dim && escape_chebysev) {
                        escape = true;
                        escape_chebysev = false;
                        escape_billiard = true;
                        stuck = false;
                    }
                    else if (stuck_try_coord == dim)
                        break;
                }
                else {
                    if (stuck && stuck_try_coord == dim) {
                        escape = true;
                        escape_chebysev = true;
                        stuck = false;
                    } else if (!stuck) {
                        stuck = true;
                        stuck_try_coord = 0;
                    }
                }
            }
            else {
                stuck = false;
                triedEscaping = false;
                escape = escape_billiard = escape_chebysev = false;
            }

            stepsSinceLastEscape++;
            step++;
        } while (step <= maxSteps || tillConvergence);
        STEPS = step - 1;

        if (verbose) std::cout << "Ended at " << step -1<< " steps " << polytope.is_in(interiorPoint) << std::endl;
        return std::pair<Point, NT>(interiorPoint, dotProduct);


    }



    /**
     * Solve the linear program
     *
     *      min cx
     *      s.t. Ax <= b
     *
     * @tparam Parameters
     * @tparam Point
     * @tparam NT
     * @param polytope
     * @param objectiveFunction
     * @param parameters
     * @param error
     * @param maxSteps
     * @param initial
     * @return
     */
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method(HPolytope<Point> polytope, VT& objectiveFunction, Parameters parameters, const NT error,
                        const unsigned int maxSteps, Point& initial) {
        normalizePolytope(polytope);

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = objectiveFunction.rows();

        SlidingWindow slidingWindow(5 + sqrt(dim));
        std::pair<Point, Point> minimizingPoints;


        // get an internal point so you can sample
        Point interiorPoint = initial;

        minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters);

        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
        slidingWindow.push(min);

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        do {
            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);

            // find where to cut the polytope
            minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters);

            min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
            slidingWindow.push(min);

            if (slidingWindow.getRelativeError() < error)
                    break;


            step++;
//            std::cout << min << ", ";


        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step<< " steps"  <<  std::endl;


        return std::pair<Point, NT>(minimizingPoints.first, objectiveFunction.dot(minimizingPoints.first.getCoefficients()));
    }


    /**
     * Solve the linear program
     *
     *      min cx
     *      s.t. Ax <= b
     *
     * @tparam Parameters
     * @tparam Point
     * @tparam NT
     * @param polytope
     * @param objectiveFunction
     * @param parameters
     * @param error
     * @param maxSteps
     * @param initial
     * @return
     */
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method_billiard(HPolytope<Point> polytope, VT& objectiveFunction, Parameters parameters, const NT error,
                                              const unsigned int maxSteps, Point& initial) {
        normalizePolytope(polytope);

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = objectiveFunction.rows();
        SlidingWindow slidingWindow(5 + sqrt(dim));
        std::pair<Point, Point> minimizingPoints;


        // get an internal point so you can sample
        Point interiorPoint = initial;

        minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters);

        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
        slidingWindow.push(min);


        double diameter;

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);
        Point previousMinimizingSecond = initial;

        do {
            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);
            diameter = euclideanDistance(previousMinimizingSecond, minimizingPoints.second);

            previousMinimizingSecond = minimizingPoints.second;

            // find where to cut the polytope
            minimizingPoints = min_rand_point_generator_billiard(polytope, objectiveFunction, interiorPoint, rnum, parameters, diameter);

            min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
            slidingWindow.push(min);

            if (slidingWindow.getRelativeError() < error)
                break;

//            std::cout << min << ", ";

            step++;


        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step<< " steps"  <<  std::endl;


        return std::pair<Point, NT>(minimizingPoints.first, objectiveFunction.dot(minimizingPoints.first.getCoefficients()));
    }

    /**
     * Solve the linear program
     *
     *      min cx
     *      s.t. Ax <= b
     *
     * with the implicit isotropization heuristic
     *
     * @tparam Parameters
     * @tparam Point
     * @tparam NT
     * @param polytope
     * @param objectiveFunction
     * @param parameters
     * @param error
     * @param maxSteps
     * @param initial
     * @return
     */
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method_sampled_covariance_matrix(HPolytope<Point> polytope, VT &objectiveFunction,
                                                             Parameters parameters, const NT error,
                                                             const unsigned int maxSteps, Point &initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m; // # of random points per phase
        unsigned int walk_len = parameters.walk_steps; // mixing time for the random walk
//        SlidingWindow slidingWindow(walk_len);
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        SlidingWindow slidingWindow(3);

        std::pair<Point, Point> minimizingPoints;

        std::list<Point> points;

        // get an internal point so you can sample
        Point interiorPoint = initial;
//        minimizingPoints = min_rand_point_generator_boundary(polytope, objectiveFunction, interiorPoint, rnum, walk_len, intersectionPoints, parameters);
        minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters, points);

        // find where to cut the polytope
        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
        slidingWindow.push(min);

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        do {
            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);

            try { // we may fail to compute square root of matrix
                MT covarianceMatrix = sampledCovarianceMatrix(points);
                points.clear();
                minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters, points, covarianceMatrix);
            }
            catch (int e) {
                // we end up here if we fail to compute the sqaure root of a matrix
                points.clear();
                minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters, points);
            }


            // check for distance between successive estimations
            min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
            slidingWindow.push(min);

            if (slidingWindow.getRelativeError() < error)
                break;

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;


        if (verbose) std::cout << "Ended at " << step  << " steps with implicit isotropization" << std::endl;
//        assert(polytope.is_in(minimizingPoints.first));

        return std::pair<Point, NT>(minimizingPoints.first, objectiveFunction.dot(minimizingPoints.first.getCoefficients()));
    }


    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method_deterministic_Chebyshev(HPolytope<Point> polytope, VT &objectiveFunction,
                                                                        Parameters parameters, const NT error,
                                                                        const unsigned int maxSteps, Point &initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m; // # of random points per phase
        unsigned int walk_len = parameters.walk_steps; // mixing time for the random walk
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        SlidingWindow slidingWindow(3);


        // get an internal point so you can sample
        Point interiorPoint = initial;

        // find where to cut the polytope
        normalizePolytope(polytope);
        escapeStep_ChebyshevCenter(polytope, interiorPoint, rnum);

        Point cutAt = interiorPoint;
        NT min = objectiveFunction.dot(cutAt.getCoefficients());
        slidingWindow.push(min);
//        std::cout << min << " " << polytope.is_in(interiorPoint) << " " << interiorPoint.getCoefficients().transpose() << "\n";

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);
        Point obj(objectiveFunction);
        Point dist(objectiveFunction*0.00001);
        escapeStep_ChebyshevCenter(polytope, interiorPoint, rnum);
        cutAt = interiorPoint + dist;

        do {
            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, cutAt);
            normalizePolytope(polytope);

            escapeStep_ChebyshevCenter(polytope, interiorPoint, rnum);
            cutAt = interiorPoint + dist;

            min = objectiveFunction.dot(cutAt.getCoefficients());

            // check for distance between successive estimations
            slidingWindow.push(min);
//            std::cout << min << " " << polytope.is_in(interiorPoint) << " " << interiorPoint.getCoefficients().transpose() << "\n";

            if (slidingWindow.getRelativeError() < error)
                break;

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;


        if (verbose) std::cout << "Ended at " << step  << " steps with deterministic Chebushev" << std::endl;
//        assert(polytope.is_in(minimizingPoints.first));

        return std::pair<Point, NT>(interiorPoint, objectiveFunction.dot(interiorPoint.getCoefficients()));
    }
}

#endif //VOLESTI_CUTTING_PLANE_H
