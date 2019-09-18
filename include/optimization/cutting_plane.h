// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Panagiotis Repouskos, as part of Google Summer of Code 2019 program.

// VolEsti is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// VolEsti is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.

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

int STEPS;

namespace optimization {


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
    MT getIsotropicQuantities(VT &c, std::list<Point> &points, std::pair<Point, Point> &pointsPair, std::list<NT> &dotProducts) {
        class std::list<Point>::iterator points_it = points.begin();
        typename std::list<NT>::iterator products_it = dotProducts.begin();

        int dim = points_it->dimension();

        NT currentMinimum = c.dot(pointsPair.second.getCoefficients());
        VT y;
        y.setZero(dim);

        for (auto p : points)
            y += p.getCoefficients() / (double) points.size();

        // compute Y
        MT Y;
        Y.setZero(dim, dim);
        VT temp(dim);

        for (auto pit = points.begin(); pit != points.end(); pit++) {
            temp = pit->getCoefficients() - y;
            Y = Y + (temp * temp.transpose());
        }

        Y = Y / (double) (points.size());

        try {
            Y = Y.sqrt();
        }
        catch (int e) {
            throw e;
        }

        return Y;
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
        assert(P.is_in(p));

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
     * @param isotropic the matrix from the implicit isotropization heuristic
     * @return (p1, p2) p1 minimizes c the most and p2 follows
     */
    template<class Polytope, class PointList, class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_iso(Polytope &P,
                                                         VT& c,
                                                         Point &p,
                                                         unsigned int rnum,
                                                         unsigned int endPointEveryKSteps,
                                                         PointList &endPoints,
                                                         Parameters &var,
                                                         MT& isotropic)
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;


        Point boundaryMin;
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

        int count = 0;

        for (unsigned int i = 1; i <= rnum; ++i) {


            if (var.cdhr_walk) {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                VT iso = isotropic.col(rand_coord);
                hit_and_run_coord_update_isotropic(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas, bpair, p1, p2, iso);
//            NT diff = bpair.first + kapa * (bpair.second - bpair.first);
//            currentProduct += c(rand_coord) * diff;
            } else
                implicit_isotropization_hit_and_run(p, P, var, p1, p2, isotropic);

            count++;
            if (count == endPointEveryKSteps) {
                endPoints.push_back(p1);
                endPoints.push_back(p2);
                count = 0;
            }

            currentProduct = p.getCoefficients().dot(c);

            getNewMinimizingPoints(p, dotProduct1, dotProduct2, min1, min2, currentProduct);
        }

        p = (min1 + min2)/2;
//        std::cout << min1.getCoefficients().dot(c) <<"\n";
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

            if (changedMin1){
                // if the new point is the new min update the boundary point

                boundaryMin2 = boundaryMin1;
                if (var.cdhr_walk) {
                    if (c(rand_coord) > 0) {
                        boundaryMin1 = p_prev;
                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.second);
                    } else {
                        boundaryMin1 = p_prev;
                        boundaryMin1.set_coord(rand_coord, boundaryMin1[rand_coord] + bpair.first);
                    }

                } else {
                    boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
                }
            }
            else if (changedMin2) {
                if (var.cdhr_walk) {
                    if (c(rand_coord) > 0) {
                        boundaryMin2 = p_prev;
                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.second);
                    } else {
                        boundaryMin2 = p_prev;
                        boundaryMin2.set_coord(rand_coord, boundaryMin2[rand_coord] + bpair.first);
                    }

                } else {
                    boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
                }
            }

//        assert(P.is_in(min1));
//        assert(P.is_in(min2));
//            std::cout << min1.getCoefficients().dot(c) << "\t" << minProduct1 <<"\n";
        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */


        // find an interior point to start the next phase
        Point _p =  min1*0.20;
        Point _p1 = min2*0.40;
        Point _p2 =  boundaryMin1*0.20;
        Point _p3 = boundaryMin2*0.20;
//        p = (min1 + min2)/2;
        p = _p + _p1 + _p2 + _p3;

//        std::cout << min1.getCoefficients().dot(c) << "\t" << P.is_in(p) << "\t" << min2.getCoefficients().dot(c) << "\n";
//        std::cout << min2.getCoefficients().dot(c) <<"\n";
//        std::cout << "b " << boundaryMin1.getCoefficients().dot(c) << " " << P.is_in(boundaryMin1) << "\n";

        return std::pair<Point, Point>(min1, min2);
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

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 1;

        SlidingWindow slidingWindow(2);
        std::pair<Point, Point> minimizingPoints;
        bool escape = false;
        int stepsSinceLastEscape = 0;

        // get an internal point so you can sample
        Point interiorPoint = initial;

        minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters);

        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
        slidingWindow.push(min);

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        do {
            std::list<Point> randPoints;

            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);

            // find where to cut the polytope
            minimizingPoints = min_rand_point_generator(polytope, objectiveFunction, interiorPoint, rnum, parameters);

            min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
            slidingWindow.push(min);

            if (slidingWindow.getRelativeError() < error) break;

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step - 1;

        if (verbose) std::cout << "Ended at " << step -1<< " steps"  <<  std::endl;

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
    std::pair<Point, NT> cutting_plane_method_isotropization(HPolytope<Point> polytope, VT &objectiveFunction,
                                                             Parameters parameters, const NT error,
                                                             const unsigned int maxSteps, Point &initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m; // # of random points per phase
        unsigned int walk_len = parameters.walk_steps; // mixing time for the random walk
//        SlidingWindow slidingWindow(walk_len);
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;

        std::pair<Point, Point> minimizingPoints;

        // the intersection points between the polytope and the lines of hit and run
        std::list<Point> intersectionPoints;

        // get an internal point so you can sample
        Point interiorPoint = initial;
        minimizingPoints = min_rand_point_generator_boundary(polytope, objectiveFunction, interiorPoint, rnum, walk_len, intersectionPoints, parameters);

        // find where to cut the polytope
        NT newMin = objectiveFunction.dot(minimizingPoints.second.getCoefficients());

        // add one more row in polytope, where we will store the current cutting plane
        // each time we cut the polytope we replace the previous cutting plane with the new one
        addRowInPolytope<Point, NT>(polytope);

        NT min;

        do {

            min = newMin;

            // cut the polytope
            cutPolytope<Point, NT>(objectiveFunction, polytope, minimizingPoints.second);

            //TODO do we need this? it works fine as is
            // we need this for arithmetic stability
//            intersectionPoints.push_back(minimizingPoints.first);
//            intersectionPoints.push_back(minimizingPoints.second);

            // sample points from polytope
            // matrix isotropic will be multiplied witch each direction vector of hit and run

            std::list<NT> intersectionPointsDotProducts;

            dotProducts(objectiveFunction, intersectionPoints, intersectionPointsDotProducts);
//            interiorPoint = getArithmeticMean(objectiveFunction, intersectionPoints, minimizingPoints.second, intersectionPointsDotProducts, sum);

            try { // we may fail to compute square root of matrix
                MT isotropic = getIsotropicQuantities<Point, NT>(objectiveFunction, intersectionPoints, minimizingPoints, intersectionPointsDotProducts);

                intersectionPoints.clear();
                minimizingPoints = min_rand_point_generator_iso(polytope, objectiveFunction, interiorPoint, rnum, walk_len, intersectionPoints, parameters, isotropic);

            }
            catch (int e) {
                // we end up here if we fail to compute the sqaure root of a matrix
                intersectionPoints.clear();
                minimizingPoints = min_rand_point_generator_boundary(polytope, objectiveFunction, interiorPoint, rnum, walk_len, intersectionPoints, parameters);
            }


            // check for distance between successive estimations
            newMin = objectiveFunction.dot(minimizingPoints.second.getCoefficients());

            if (relative_error(min, objectiveFunction.dot(minimizingPoints.first.getCoefficients())) < error) break;


            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;


        if (verbose) std::cout << "Ended at " << step  << " steps with implicit isotropization" << std::endl;
//        assert(polytope.is_in(minimizingPoints.first));

        return std::pair<Point, NT>(minimizingPoints.first, objectiveFunction.dot(minimizingPoints.first.getCoefficients()));
    }
}

#endif //VOLESTI_CUTTING_PLANE_H
