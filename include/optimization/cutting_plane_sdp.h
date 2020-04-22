//
// Created by panagiotis on 24/5/2019.
//

#ifndef VOLESTI_CUTTING_PLANE_SDP_H
#define VOLESTI_CUTTING_PLANE_SDP_H


#include "spectrahedron.h"
#include "Eigen"
#include <list>
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
#include "samplers.h"
#include "interior_point.h"
#include <queue>
#include "heuristics.h"

//int STEPS;
extern int BOUNDARY_CALLS;
extern double ORACLE_TIME;
extern double REFLECTION_TIME;
namespace optimization {


//    const double ZERO = 0.000000000001;

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
//    double relative_error(double approx, double exact) {
//        return abs((exact - approx) / exact);
//    }


    /**
     * A class to save successive approximation values
     */
//    class SlidingWindow {
//    public:
//        std::queue<double> evals;
//        int evalsSize; // how many values to store
//        int count;
//        double min;
//
//        SlidingWindow(int evalSize) {
//            this->evalsSize = evalSize;
//            count = 0;
//        }
//
//        void push(double eval) {
//            if (count >= evalsSize) {
//                evals.pop();
//            } else
//                count++;
//
//            evals.push(eval);
//            min = eval;
//        }
//
//        double getRelativeError() {
//            if (count < evalsSize)
//                return 1;
//
//            return abs((min - evals.front()) / min);
//        }
//
//        void half() {
//            for (int i = 0; i < evals.size() / 2; i++) {
//                evals.pop();
//                count--;
//            }
//        }
//    };


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
//    template<class Point, typename NT>
//    void
//    getNewMinimizingPoints(const Point &p, NT &dotProduct1, NT &dotProduct2, Point &min1, Point &min2, NT newProduct,
//                           bool &changedMin1, bool &changedMin2) {
//
//        if (newProduct < dotProduct2) {
//            if (newProduct < dotProduct1) {
//                dotProduct2 = dotProduct1;
//                min2 = min1;
//                dotProduct1 = newProduct;
//                min1 = p;
//                changedMin1 = true;
//            } else {
//                dotProduct2 = newProduct;
//                min2 = p;
//                changedMin2 = true;
//            }
//        }
//
//    }

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
//    template<class Point, typename NT>
//    void
//    getNewMinimizingPoints(const Point &p, NT &dotProduct1, NT &dotProduct2, Point &min1, Point &min2, NT newProduct) {
//
//        if (newProduct < dotProduct2) {
//            if (newProduct < dotProduct1) {
//                dotProduct2 = dotProduct1;
//                min2 = min1;
//                dotProduct1 = newProduct;
//                min1 = p;
//            } else {
//                dotProduct2 = newProduct;
//                min2 = p;
//            }
//        }
//
//    }

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
    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_spectrahedron(Spectrahedron &spectrahedron,
                                                     VT &c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                     std::list<Point>& points)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);

        hit_and_run(p, spectrahedron, var);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        hit_and_run(p, spectrahedron, var);

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

        for (unsigned int i = 1; i <= rnum; ++i) {

            hit_and_run(p, spectrahedron, var, p1, p2);
            newProduct = p.getCoefficients().dot(c);
            points.push_back(p);


            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

            if (changedMin1) {
                // if the new point is the new min update the boundary point

                boundaryMin2 = boundaryMin1;
                boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;

            } else if (changedMin2) {
                boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;

            }

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

//TODO done need boundary perhaps?
        // find an interior point to start the next phase
        Point _p = min1 * 0.50;
        Point _p1 = min2 * 0.50;
        p = _p + _p1;

        return std::pair<Point, Point>(min1, min2);
    }


    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_spectrahedron(Spectrahedron &spectrahedron,
                                                     VT &c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                     VT &a,
                                                     double b,
                                                     std::list<Point>& points)
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);

        hit_and_run(p, spectrahedron, var, a, b);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        hit_and_run(p, spectrahedron, var, a, b);

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
        int addPointEverySteps = rnum / (1000 + dim*sqrt(dim));

        if (addPointEverySteps == 0)
            addPointEverySteps = 1;

        for (unsigned int i = 1; i <= rnum; ++i) {

            // get next point

            hit_and_run(p, spectrahedron, var, p1, p2, a, b);
            newProduct = p.getCoefficients().dot(c);

            if (i % addPointEverySteps == 0)
                points.push_back(p);

            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

        // find an interior point to start the next phase
        Point _p = min1 * 0.50;
        Point _p1 = min2 * 0.50;
        p = _p + _p1;

        return std::pair<Point, Point>(min1, min2);
    }

    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_spectrahedron(Spectrahedron &spectrahedron,
                                                     VT &c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                     VT &a,
                                                     double b,
                                                     std::list<Point>& points,
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

        hit_and_run_sampled_covariance_matrix(p, spectrahedron, var, a, b, covarianceMatrix);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        hit_and_run_sampled_covariance_matrix(p, spectrahedron, var, a, b, covarianceMatrix);


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

        int addPointEverySteps = rnum / (1000 + dim*sqrt(dim));

        if (addPointEverySteps == 0)
            addPointEverySteps = 1;

        for (unsigned int i = 1; i <= rnum; ++i) {

            // get next point

            hit_and_run_sampled_covariance_matrix(p, spectrahedron, var, p1, p2, a, b, covarianceMatrix);
            newProduct = p.getCoefficients().dot(c);

            if (i % addPointEverySteps == 0)
                points.push_back(p);


            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);


        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

//TODO done need boundary perhaps?
        // find an interior point to start the next phase
        Point _p = min1 * 0.50;
        Point _p1 = min2 * 0.50;
        p = _p + _p1;// + _p2 + _p3;

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
    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_spectrahedron(Spectrahedron &spectrahedron,
                                                     const VT &c,
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

        hit_and_run(p, spectrahedron, var);


        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        hit_and_run(p, spectrahedron, var);

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

        for (unsigned int i = 1; i <= rnum; ++i) {

            hit_and_run(p, spectrahedron, var, p1, p2);
            newProduct = p.getCoefficients().dot(c);



            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

//            if (changedMin1) {
                // if the new point is the new min update the boundary point
//
//                boundaryMin2 = boundaryMin1;
//                boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//
//            } else if (changedMin2) {
//                boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
//
//            }

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

//TODO done need boundary perhaps?
        // find an interior point to start the next phase
        Point _p = min1 * 0.50;
        Point _p1 = min2 * 0.50;
        p = _p + _p1;

        return std::pair<Point, Point>(min1, min2);
    }


    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_spectrahedron(Spectrahedron &spectrahedron,
                                                     VT &c,
                                                     Point &p,   // a point to start
                                                     unsigned int rnum,
                                                     Parameters &var,
                                                     VT &a,
                                                     double b)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);

        Point p1(dim), p2(dim), min1(dim), min2(dim);

        hit_and_run(p, spectrahedron, var, a, b);


        // get the first two points
        min1 = p;
        NT minProduct1 = min1.getCoefficients().dot(c);

        hit_and_run(p, spectrahedron, var, a, b);

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

        for (unsigned int i = 1; i <= rnum; ++i) {

            // get next point

            hit_and_run(p, spectrahedron, var, p1, p2, a, b);
            newProduct = p.getCoefficients().dot(c);



            // get new minimizing point
            bool changedMin1 = false;
            bool changedMin2 = false;

            getNewMinimizingPoints(p, minProduct1, minProduct2, min1, min2, newProduct, changedMin1, changedMin2);

            if (changedMin1) {
                // if the new point is the new min update the boundary point

                boundaryMin2 = boundaryMin1;
                boundaryMin1 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;
            } else if (changedMin2) {
                boundaryMin2 = p1.getCoefficients().dot(c) < p2.getCoefficients().dot(c) ? p1 : p2;

            }
        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

//TODO done need boundary perhaps?
        // find an interior point to start the next phase
        Point _p = min1 * 0.50;
        Point _p1 = min2 * 0.50;
        p = _p + _p1;// + _p2 + _p3;

        return std::pair<Point, Point>(min1, min2);
    }



    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    cutting_plane_method_sampled_covariance_matrix(Spectrahedron &spectrahedron, VT &objectiveFunction, Parameters parameters, const NT error,
                         const unsigned int maxSteps, Point &initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = objectiveFunction.rows();

        SlidingWindow slidingWindow(5 + sqrt(dim));
        std::pair<Point, Point> minimizingPoints;
        std::list<Point> points;

        // get an internal point so you can sample
        Point interiorPoint = initial;


        minimizingPoints = min_rand_point_generator_spectrahedron(spectrahedron, objectiveFunction, interiorPoint, rnum, parameters, points);

        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
        slidingWindow.push(min);

        VT a  = objectiveFunction;
        double b = objectiveFunction.dot(minimizingPoints.second.getCoefficients());

        do {
            std::list<Point> randPoints;

            try {
                MT covarianceMatrix = sampledCovarianceMatrix(points);
                points.clear();
//                std::cout << covarianceMatrix << "\n";
                minimizingPoints = min_rand_point_generator_spectrahedron(spectrahedron, objectiveFunction, interiorPoint, rnum,
                                                            parameters, a, b, points, covarianceMatrix);
            }
            catch (int e) {
//                if (verbose) std::cout << "Failed to compute covariance matrix - step " << step << "\n";
                points.clear();
                minimizingPoints = min_rand_point_generator_spectrahedron(spectrahedron, objectiveFunction, interiorPoint, rnum,
                                                            parameters, a, b, points);
            }


            min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
            b = min;
            slidingWindow.push(min);

            VT coeffs = interiorPoint.getCoefficients();
            if (slidingWindow.getRelativeError() < error || spectrahedron.isSingular(coeffs))
                    break;

            step++;

//            LMI lmi;
//            lmi = spectrahedron.getLMI();
//            VT coeffs = minimizingPoints.first.getCoefficients();
//            std::cout << min << ", ";
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;


        return std::pair<Point, NT>(minimizingPoints.first,
                                    objectiveFunction.dot(minimizingPoints.first.getCoefficients()));
    }

    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    cutting_plane_method(Spectrahedron &spectrahedron, VT &objectiveFunction, Parameters parameters, const NT error,
                         const unsigned int maxSteps, Point &initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = objectiveFunction.rows();

        SlidingWindow slidingWindow(5 + sqrt(dim));
        std::pair<Point, Point> minimizingPoints;


        // get an internal point so you can sample
        Point interiorPoint = initial;

        minimizingPoints = min_rand_point_generator_spectrahedron(spectrahedron, objectiveFunction, interiorPoint, rnum, parameters);

        NT min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
        slidingWindow.push(min);

        VT a  = objectiveFunction;
        double b = objectiveFunction.dot(minimizingPoints.second.getCoefficients());

        do {
            std::list<Point> randPoints;

            // find where to cut the polytope
            minimizingPoints = min_rand_point_generator_spectrahedron(spectrahedron, objectiveFunction, interiorPoint, rnum,
                                                        parameters, a, b);

            min = objectiveFunction.dot(minimizingPoints.second.getCoefficients());
            b = min;
            slidingWindow.push(min);

            VT coeffs = interiorPoint.getCoefficients();
            if (slidingWindow.getRelativeError() < error || spectrahedron.isSingular(coeffs))
                break;


            step++;

//            std::cout << min << " " << step << "+++++++++++++++++++++++++++++++++++++++++++\n";
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;


        return std::pair<Point, NT>(minimizingPoints.first,
                                    objectiveFunction.dot(minimizingPoints.first.getCoefficients()));
    }

    template <class Point, typename NT>
    void checkNewMin(const NT& newVal,const Point& newPoint, NT& min1, Point& minP1, NT& min2, Point& minP2) {
        if (newVal < min1) {
            min2 = min1;
            min1 = newVal;
            minP2 = minP1;
            minP1 = newPoint;
        }
        else if (newVal < min2) {
            min2 = newVal;
            minP2 = newPoint;
        }
    }

    template <class Point, typename NT>
    void findMinInList(Point& retMin, const Point& objectiveFunction, const std::list<Point>& pointlist) {
        NT min = pointlist.front().dot(objectiveFunction);
        retMin = pointlist.front();

        for (auto point : pointlist) {
            NT newProduct = point.dot(objectiveFunction);
            if (newProduct < min) {
                min = newProduct;
                retMin = point;
            }
        }
    }

    template <class Point, class Pointlist, typename NT>
    void checkNewMinInList(Point& minP1, Point& minP2, NT& min1, NT& min2, const Point& objectiveFunction, const Pointlist& intermediateBilliardPoints) {
        if (intermediateBilliardPoints.empty())
            return;

        Point minInList;
        findMinInList<Point, NT>(minInList, objectiveFunction, intermediateBilliardPoints);
        checkNewMin(minInList.dot(objectiveFunction), minInList, min1, minP1, min2, minP2);
    }


    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_cdhr(Spectrahedron &spectrahedron,
                                                              const Point& objectiveFunction,
                                                              Point &p,   // a point to start
                                                              const unsigned int& rnum,
                                                              const Point& a,
                                                              const double& b,
                                                              Parameters &var,
                                                              Spectrahedron::BoundaryOracleCDHRSettings<Point>& settings)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        RNGType &rng = var.rng;
        boost::random::uniform_int_distribution<> uidist(0, spectrahedron.getLMI().getDim() - 1);

        Point minP1, minP2;
        double lambda;

        hit_and_run_coord_update(p, spectrahedron, uidist(rng), a, b, var, settings, lambda);
        settings.first = false;
        minP1 = p;
        hit_and_run_coord_update(p, spectrahedron, uidist(rng), a, b, var, settings, lambda);
        minP2 = p;

        NT min1 = minP1.dot(objectiveFunction);
        NT min2 = minP2.dot(objectiveFunction);

        if (min1 > min2) {
            std::swap(min1, min2);
            std::swap(minP1, minP2);
        }

        double currentProduct = p.dot(objectiveFunction);

        for (unsigned int i = 3; i <= rnum ; ++i) {

            int coord = uidist(rng);
            hit_and_run_coord_update(p, spectrahedron, coord, a, b, var, settings, lambda);
            currentProduct += objectiveFunction[coord] * lambda;

            checkNewMin(currentProduct, p, min1, minP1, min2, minP2);

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

        p = 0.5*minP1 + 0.5*minP2;// + _p2 + _p3;
        settings.first = true;

        return std::pair<Point, Point>(minP1, minP2);
    }


    template<class Parameters, class Point>
    std::pair<Point, Point> min_rand_point_generator_billiard(Spectrahedron &spectrahedron,
                                                              Point& objectiveFunction,
                                                              Point &p,   // a point to start
                                                              unsigned int rnum,
                                                              const Point& a,
                                                              const double& b,
                                                              Parameters &var,
                                                              double diameter,
                                                              Spectrahedron::BoundaryOracleBilliardSettings<Point>& settings)
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

//        std::list<Point> intermediateBilliardPoints;

        Point minP1, minP2;

        if (rnum == 1) {
            while (billiard_walk(spectrahedron, p, 1.5*diameter,  var, a, b, settings, true) < 1);
//            settings.first = false;
//            double lambda = 0.005;
//            std::pair<double, bool> pair;
//            pair = spectrahedron.boundaryOracleBilliard(p.getCoefficients(), objectiveFunction.getCoefficients(), a.getCoefficients(), b, settings);
//            minP1 = p;
//            minP2 = p + pair.first*lambda*objectiveFunction;
//            p = 0.5*minP1 + 0.5*minP2;// + _p2 + _p3;
//            settings.first = true;
            return std::pair<Point, Point>(p, p);
        }

        billiard_walk(spectrahedron, p, diameter,  var, a, b, settings);
        settings.first = false;

        minP1 = p;
        billiard_walk(spectrahedron, p, diameter,  var, a, b, settings);
        minP2 = p;

        NT min1 = minP1.dot(objectiveFunction);
        NT min2 = minP2.dot(objectiveFunction);

        if (min1 > min2) {
            std::swap(min1, min2);
            std::swap(minP1, minP2);
        }

//        checkNewMinInList(minP1,minP2, min1, min2, objectiveFunction, intermediateBilliardPoints);
//        intermediateBilliardPoints.clear();


        for (unsigned int i = 3; i <= rnum ; ++i) {

            billiard_walk(spectrahedron, p, diameter,  var, a, b,  settings);
//            std::cout << "size " << intermediateBilliardPoints.size() << "\n";
            NT newProduct = p.dot(objectiveFunction);
            checkNewMin(newProduct, p, min1, minP1, min2, minP2);

//            checkNewMinInList(minP1,minP2, min1, min2, objectiveFunction, intermediateBilliardPoints);
//            intermediateBilliardPoints.clear();

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

        // find an interior point to start the next phase
//        Point _p =  0.5*minP1;
//        Point _p1 = 0.5*minP2;

        p = 0.5*minP1 + 0.5*minP2;// + _p2 + _p3;
        settings.first = true;
//        std::cout << "res " << minP1.dot(objectiveFunction) << " " << minP2.dot(objectiveFunction) << " " << euclideanDistance(minP1, p) << "\n";
        return std::pair<Point, Point>(minP1, minP2);
    }

    template <class Point>
    double getDiameterLength(const Point& p1, const Point& p2, double maxSegmentLength, double diameter) {
        double distance = euclideanDistance(p1, p2);
        double retVal = maxSegmentLength > distance ? maxSegmentLength : distance;
        if (retVal < diameter)
            retVal = diameter;

        return (1.5 + (double) pow((double)p1.dimension(), 0.25) / 5) * retVal;
    }

    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    cutting_plane_method_billiard(Spectrahedron &spectrahedron, Point &objectiveFunction, Parameters parameters, const NT error,
                         const unsigned int maxSteps, Point &initial) {

        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = objectiveFunction.dimension();
        BOUNDARY_CALLS = ORACLE_TIME = REFLECTION_TIME = 0;

        SlidingWindow slidingWindow(6 + sqrt(dim));
        SlidingWindow changeWindow(3 + 2*sqrt(dim)/3);

        std::pair<Point, Point> minimizingPoints;
        Spectrahedron::BoundaryOracleBilliardSettings<Point> billiardSettings(spectrahedron.getLMI().getMatricesDim());
        Spectrahedron::BoundaryOracleCDHRSettings<Point> cdhrSettings(spectrahedron.getLMI().getMatricesDim());

        //prepare cutting plane
        Point a = objectiveFunction;
        double a_norm = a.normalize();
        double b;
        double diameter;

        // get an internal point so you can sample
        Point interiorPoint = initial;
        minimizingPoints = min_rand_point_generator_cdhr(spectrahedron, objectiveFunction, interiorPoint, 2*dim + 10, a, b, parameters, cdhrSettings);

        NT min = objectiveFunction.dot(minimizingPoints.second);
        slidingWindow.push(min);
        changeWindow.push(min);

        b = objectiveFunction.dot(minimizingPoints.second) / a_norm;
        cdhrSettings.hasCuttingPlane = true;


        Point previousMinimizingSecond;
        diameter = cdhrSettings.maxSegmentLength();
        bool useBilliard = false;
        int previous = 0;
        double lastPrevious = min;

        do {
            previousMinimizingSecond = minimizingPoints.second;

            // find where to cut the polytope

            cdhrSettings.resetMaxSegmentLength();

//            if (changeWindow.getRelativeError() > 0.0001) {
            if (!useBilliard) {
                minimizingPoints = min_rand_point_generator_cdhr(spectrahedron, objectiveFunction, interiorPoint, dim*1.5, a, b, parameters, cdhrSettings);
            }
            else {
                billiardSettings.resetMaxSegmentLength();
                minimizingPoints = min_rand_point_generator_billiard(spectrahedron, objectiveFunction, interiorPoint, rnum, a, b, parameters, diameter, billiardSettings);
            }

            min = objectiveFunction.dot(minimizingPoints.second);
            b = min/a_norm;
//            std::cout << BOUNDARY_CALLS - previous << " " << min << "\n";
//            previous = BOUNDARY_CALLS;



            if (!useBilliard) {
                if (changeWindow.getRelativeErrorLastRound(min) < 0.001) {
//                    std::cout << changeWindow.getRelativeErrorLastRound(min) << " " << step
//                              << " -------------------------------------\n";
                    useBilliard = true;
                    step++;
                    changeWindow.push(min);
                    slidingWindow.push(min);
                    lastPrevious = min;
                    continue;
                }
            }

            changeWindow.push(min);
            slidingWindow.push(min);

            if (slidingWindow.getRelativeError() < error)
                break;

            if (useBilliard)
                diameter = getDiameterLength(previousMinimizingSecond, minimizingPoints.first, billiardSettings.maxSegmentLength(), cdhrSettings.maxSegmentLength());

            step++;
//            std::cout << changeWindow.getRelativeError() << "\n";
            if ((changeWindow.getRelativeError() > 0.01 || relative_error(lastPrevious, min) > 0.01) && useBilliard) {
//                std::cout << changeWindow.getRelativeErrorLastRound(min) << " ++++++++++++++++++++++++++++\n";
                useBilliard = false;
            }

        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;


        return std::pair<Point, NT>(minimizingPoints.first,
                                    objectiveFunction.dot(minimizingPoints.first));
    }

}


#endif //VOLESTI_CUTTING_PLANE_SDP_H
