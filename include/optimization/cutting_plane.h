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

namespace optimization {



    typedef double NT_MATRIX;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT_MATRIX, Eigen::Dynamic, 1> VT;



    /**
     * Compute the cutting plane c (x - point) <= 0  ==>  cx <= c point
     * and add it as a new constraint in the existing polytope, in place of its last one, which will be redundant.
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

        //add cx in last row of A
        long j, i;

        for (j = 0, i = polytope.num_of_hyperplanes() - 1 ; j < dim; j++) {
            polytope.put_mat_coeff(i, j, c(j));
        }

        //add  < c,  point >  in last row of b
        NT _b = c.dot(point.getCoefficients());
        polytope.put_vec_coeff(polytope.num_of_hyperplanes() - 1, _b);
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
     * Computes the quantites y = (1/N) * Sum {p | p in points}
     * and Y = (1/N) * Sum [(p-y)* (p -y)^T], p in points
     *
     * These quantities will be used to set the direction vector of hit and run
     *
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param points a collection of points
     * @param dotProducts a list with the dot products of c with all elements of points
     * @param sum of elems in points
     * @return sqrt(Y)
     */
    template<class Point, typename NT>
    MT getIsotropicQuantities(VT &c, std::list<Point> &points, std::pair<Point, Point> &pointsPair, std::list<NT> &dotProducts, VT &sum) {
        class std::list<Point>::iterator points_it = points.begin();
        typename std::list<NT>::iterator products_it = dotProducts.begin();

        int dim = points_it->dimension();

        NT currentMinimum = c.dot(pointsPair.second.getCoefficients());
        VT y = sum / (double) points.size();

        // compute Y
        MT Y;
        Y.setZero(dim, dim);
        VT temp(dim);

        for (auto pit = points.begin(); pit != points.end(); pit++) {
            temp = pit->getCoefficients() - y;
            Y = Y + (temp * temp.transpose());
        }

        Y = Y / (double) points.size();

        try {
            Y = Y.sqrt();
        }
        catch (int e) {
            throw e;
        }

        return Y;
    }


    /**
     * Solve the linear program
     *
     *      min cx
     *      s.t. Ax <= b
     *
     * using the optimization algorithm in http://www.optimization-online.org/DB_FILE/2008/12/2161.pdf
     *
     * @tparam Parameters struct vars
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param lp
     * @param parameters An instance of struct vars
     * @param error how much distance between two successive estimations before we stop
     * @param masSteps maximum number of steps
     * @return A pair of the point that minimizes the object function and the minimum value
     */
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method(HPolytope<Point> polytope, VT& objectiveFunction, Parameters parameters, const NT error,
                        const unsigned int maxSteps, Point& initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        unsigned int walk_len = parameters.walk_steps;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 1;


        std::pair<Point, Point> minimizingPoints;
        std::list<Point> randPoints;
        // the intersection points between the polytope and the lines of hit and run
        std::list<Point> intersectionPoints;

        // get an internal point so you can sample
        auto t1 = std::chrono::steady_clock::now();

        // get an internal point so you can sample
        Point interiorPoint = initial;

        rand_point_generator(polytope, interiorPoint, rnum, walk_len, randPoints, intersectionPoints, parameters);

        // find where to cut the polytope
        minimizingPoints = getPairMinimizingPoint<Point, NT>(objectiveFunction, randPoints);
        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        do {
            std::list<Point> randPoints;

            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);

            // we need this for arithmetic stability
            intersectionPoints.push_back(minimizingPoints.first);
            intersectionPoints.push_back(minimizingPoints.second);

            // sample points from polytope

            // use this point as the next starting point for sampling
            interiorPoint = getArithmeticMean<Point, NT>(objectiveFunction, intersectionPoints, minimizingPoints.second);

            intersectionPoints.clear();
            rand_point_generator(polytope, interiorPoint, rnum, walk_len, randPoints, intersectionPoints, parameters);

            // find where to cut the polytope
            minimizingPoints = getPairMinimizingPoint<Point, NT>(objectiveFunction, randPoints);

            NT newMin = objectiveFunction.dot(minimizingPoints.first.getCoefficients());
            NT distance = abs(newMin - min);
            min = newMin;

            if (distance < error) break;

            step++;
        } while (step <= maxSteps || tillConvergence);


        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;

        return std::pair<Point, NT>(minimizingPoints.first, min);
    }


    /**
     * Solve the linear program
     *
     *      min cx
     *      s.t. Ax <= b
     *
     * using the optimization algorithm in http://www.optimization-online.org/DB_FILE/2008/12/2161.pdf
     *
     * @tparam Parameters struct vars
     * @tparam Point class Point
     * @tparam NT The numeric type
     * @param lp
     * @param parameters An instance of struct vars
     * @param error how much distance between two successive estimations before we stop
     * @param masSteps maximum number of steps
     * @return A pair of the point that minimizes the object function and the minimum value
     */
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT> cutting_plane_method_isotropic(HPolytope<Point> polytope, VT& objectiveFunction, Parameters parameters, const NT error,
                        const unsigned int maxSteps, Point& initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        unsigned int walk_len = parameters.walk_steps;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 1;

        std::pair<Point, Point> minimizingPoints;
        std::list<Point> randPoints;
        // the intersection points between the polytope and the lines of hit and run
        std::list<Point> intersectionPoints;

        // get an internal point so you can sample
        Point interiorPoint = initial;
        rand_point_generator(polytope, interiorPoint, rnum, walk_len, randPoints, intersectionPoints, parameters);

        // find where to cut the polytope
        minimizingPoints = getPairMinimizingPoint<Point, NT>(objectiveFunction, randPoints);
        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        do {
            std::list<Point> randPoints;

            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);

            // we need this for arithmetic stability
            intersectionPoints.push_back(minimizingPoints.first);
            intersectionPoints.push_back(minimizingPoints.second);

            // sample points from polytope
            // matrix isotropic will be multiplied witch each direction vector of hit and run

            std::list<NT> intersectionPointsDotProducts;
            VT sum;

            dotProducts(objectiveFunction, intersectionPoints, intersectionPointsDotProducts);
            interiorPoint = getArithmeticMean(objectiveFunction, intersectionPoints, minimizingPoints.second, intersectionPointsDotProducts, sum);

            try { // we may fail to compute square root of matrix
                MT isotropic = getIsotropicQuantities<Point, NT>(objectiveFunction, intersectionPoints, minimizingPoints,
                                                                 intersectionPointsDotProducts, sum);

                intersectionPoints.clear();
                smart_rand_point_generator(polytope, interiorPoint, rnum, walk_len, randPoints,
                                           intersectionPoints, parameters, isotropic);
            }
            catch (int e) {
                intersectionPoints.clear();
                rand_point_generator(polytope, interiorPoint, rnum, walk_len, randPoints, intersectionPoints,
                                     parameters);
            }


            // find where to cut the polytope
            minimizingPoints = getPairMinimizingPoint<Point, NT>(objectiveFunction, randPoints);

            // check for distance between successive estimations
            NT newMin = objectiveFunction.dot(minimizingPoints.first.getCoefficients());
            NT distance = abs(newMin - min);
            min = newMin;

            if (distance < error) break;

            step++;
        } while (step <= maxSteps || tillConvergence);


        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;

        return std::pair<Point, NT>(minimizingPoints.first, min);
    }
}

#endif //VOLESTI_CUTTING_PLANE_H
