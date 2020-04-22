//
// Created by panagiotis on 14/8/2019.
//

#ifndef VOLESTI_SIMULATED_ANNEALING_SDP_H
#define VOLESTI_SIMULATED_ANNEALING_SDP_H

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "truncated_exponential.h"
#include "spectrahedron.h"
#include "samplers.h"
#include "Eigen"

namespace optimization {

    template<class Parameters, class Point>
    double maxDistanceFromPoint(Spectrahedron &spectrahedron,
                                                     Point &p,   // a point to start
                                                     Parameters &var)  // constants for volume
    {

        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;

        int dim = p.dimension();

        // init the walks
        RNGType &rng = var.rng;

        Point p1 = p;
        hit_and_run(p, spectrahedron, var);
        double maxDist = euclideanDistance(p1, p);

        for (unsigned int i = 1; i <= 10000 + 2*dim; ++i) {

            hit_and_run(p, spectrahedron, var);

            double newDist = euclideanDistance(p, p1);

            if (newDist > maxDist)
                maxDist = newDist;

        } /*  for (unsigned int i = 1; i <= rnum ; ++i)  */

        return maxDist;
    }

    template <class Point, class Parameters>
    void min_hit_and_run_Boltzmann(Point &p,
                                   Spectrahedron &spectrahedron,
                                   Parameters &var,
                                   Point& BoltzmannDirection,
                                   double BoltzmannParameter,
                                   const MT& choleskyDecomp,
                                   Point &minPoint,
                                   double& minValue,
                                   double& avgMinPerPhase) {
        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;
        unsigned int n =p.dimension();
        RNGType &rng = var.rng;

        Point l = get_direction<RNGType, Point, NT>(n, choleskyDecomp);

        VT pVT = p.getCoefficients();
        VT lVT = l.getCoefficients();
        std::pair <NT, NT> dbpair = spectrahedron.boundaryOracle(pVT, lVT);
        NT min_plus = dbpair.first;
        NT max_minus = dbpair.second;
        Point b1 = (min_plus * l) + p;
        Point b2 = (max_minus * l) + p;
        double c1 = BoltzmannDirection.dot(b1);
        double c2 = BoltzmannDirection.dot(b2);

        double lambda;

        if (c1 > c2) {
            lambda = texp((c1 - c2) / BoltzmannParameter, 0, min_plus - max_minus, rng);
            p = b2;
            avgMinPerPhase += c2;


            if (c2 < minValue) {
                minPoint = b2;
                minValue = c2;
            }
        }
        else {
            lambda = -texp((c2 - c1) / BoltzmannParameter, 0, min_plus - max_minus, rng);
            p = b1;
            avgMinPerPhase += c1;

            if (c1 < minValue) {
                minPoint = b1;
                minValue = c1;
            }
        }

        p = (lambda * l) + p;
    }

    template<class Parameters, class Point>
    void min_rand_point_generator_Boltzmann_window(Spectrahedron& spectrahedron,
                                                   Point& c,
                                                   Point &p,   // a point to start
                                                   double error,
                                                   Parameters &var,
                                                   double temperature,
                                                   const MT& covariance_matrix,
                                                   Point &minPoint,
                                                   double& minValue,
                                                   std::list<Point>& points,
                                                   double& avgMinPerPhase) {


        int d = p.dimension();
        double sqrtd = sqrt(d);
        int pointsSize = 1000 + 10*d*sqrtd;
        int slidingWindowSize = 1000 + sqrtd*d*d;
        int pickEverySteps = slidingWindowSize / pointsSize;

        if (pickEverySteps == 0)
            pickEverySteps = 1;

        SlidingWindow slidingWindow(slidingWindowSize);
        int count = 1;

        min_hit_and_run_Boltzmann(p, spectrahedron, var, c, temperature, covariance_matrix, minPoint, minValue, avgMinPerPhase);
        slidingWindow.push(avgMinPerPhase);
        points.push_back(p);

        while (slidingWindow.getRelativeError() > 0.00001) {
            min_hit_and_run_Boltzmann(p, spectrahedron, var, c, temperature, covariance_matrix, minPoint, minValue, avgMinPerPhase);
            count++;
            slidingWindow.push(avgMinPerPhase / (double) count);
            if (count % pickEverySteps == 0) points.push_back(p);
        }

//        std::cout << "-"<<count ;
    }

    template <class Point, class Parameters>
    void min_hit_and_run_Boltzmann(Point &p,
                                   Spectrahedron &spectrahedron,
                                   Parameters &var,
                                   Point& BoltzmannDirection,
                                   double BoltzmannParameter,
                                   Point &minPoint,
                                   double& minValue) {
        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;
        unsigned int n =p.dimension();
        RNGType &rng = var.rng;

        Point l = get_direction<RNGType, Point, NT>(n);
        VT pVT = p.getCoefficients();
        VT lVT = l.getCoefficients();
        std::pair <NT, NT> dbpair = spectrahedron.boundaryOracle(pVT, lVT);
        NT min_plus = dbpair.first;
        NT max_minus = dbpair.second;
        Point b1 = (min_plus * l) + p;
        Point b2 = (max_minus * l) + p;
        double c1 = BoltzmannDirection.dot(b1);
        double c2 = BoltzmannDirection.dot(b2);

        double lambda;

        if (c1 > c2) {
            lambda = texp((c1 - c2) / BoltzmannParameter, 0, min_plus - max_minus, rng);
            p = b2;

            if (c2 < minValue) {
                minPoint = b2;
                minValue = c2;
            }
        }
        else {
            lambda = -texp((c2 - c1) / BoltzmannParameter, 0, min_plus - max_minus, rng);
            p = b1;

            if (c1 < minValue) {
                minPoint = b1;
                minValue = c1;
            }
        }

        p = (lambda * l) + p;
    }

    template<class Parameters, class Point>
    void min_rand_point_generator_Boltzmann(Spectrahedron &spectrahedron,
                                            Point& c,
                                            Point &p,   // a point to start
                                            unsigned int walk_length,
                                            Parameters &var,
                                            double temperature,
                                            Point &minPoint,
                                            double& minValue,
                                            std::list<Point>& points) {

        // begin sampling
        for (unsigned int i = 0 ; i < walk_length ; ++i) {
            min_hit_and_run_Boltzmann(p, spectrahedron, var, c, temperature, minPoint, minValue);
            points.push_back(p);
        }
    }

    template<class Parameters, class Point>
    void min_rand_point_generator_Boltzmann(Spectrahedron &spectrahedron, Point& objectiveFunction, Point &p, const unsigned int& walk_length, double che_rad,
                                            Parameters &var,  double& temperature,  Point &minPoint, double& minValue, Spectrahedron::BoundaryOracleBoltzmannHMCSettings<Point>& settings) {

        HMC_boltzmann_reflections(spectrahedron, p, che_rad, var, objectiveFunction, temperature, settings);
        minPoint = p;
        minValue = p.dot(objectiveFunction);

        // begin sampling
        for (unsigned int i = 0 ; i < walk_length ; ++i) {
            HMC_boltzmann_reflections(spectrahedron, p, che_rad, var, objectiveFunction, temperature, settings);

            double  temp = p.dot(objectiveFunction);

            if (temp < minValue) {
                minValue = temp;
                minPoint = p;
            }
        }
    }


    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    simulated_annealing_efficient_covariance(Spectrahedron &spectrahedron, Point &objectiveFunction,
                                             Parameters parameters, const NT error,
                                             const unsigned int maxSteps, Point &initial) {


        VT obj = objectiveFunction.getCoefficients();
        double normal = obj.norm();
        obj.normalize();
        Point objFunction = Point(obj);
        bool verbose = parameters.verbose;
        unsigned int rnum = parameters.m;
        unsigned int walk_length = parameters.walk_steps;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = initial.dimension();
        std::list<Point> points;

        SlidingWindow slidingWindowComputeCovariance(3);
        SlidingWindow slidingWindowStop(5 + sqrt(dim));

        Point interiorPoint = initial;
        MT covarianceMatrix;

        double tempDescentFactor = 1 - 1 / (double) std::sqrt(dim);
        double temperature = 5;//TODO
        double temperature_threshold = 0.000001 / dim;
        double avgMinPerPhase;

        double min = interiorPoint.dot(objFunction);
        Point minPoint = interiorPoint;

        min_rand_point_generator_Boltzmann(spectrahedron, objFunction, interiorPoint, walk_length, parameters, temperature,
                                           minPoint, min, points);

        compute_covariance_matrix(points, covarianceMatrix);
        choleskyDecomposition(covarianceMatrix);


        do {

            if (temperature > temperature_threshold) {
                temperature *= tempDescentFactor;
//            else temperature *= 100;
                points.clear(); //TODO this may become too big
            }

            avgMinPerPhase = 0;
//            min_rand_point_generator_Boltzmann(polytope, objFunction, interiorPoint, walk_length, parameters, temperature, covarianceMatrix, minPoint, min, points, avgMinPerPhase);
            min_rand_point_generator_Boltzmann_window(spectrahedron, objFunction, interiorPoint, error, parameters,
                                                      temperature, covarianceMatrix, minPoint, min, points,
                                                      avgMinPerPhase);

            slidingWindowComputeCovariance.push(min);
            slidingWindowStop.push(avgMinPerPhase / (double) walk_length);


            if (slidingWindowComputeCovariance.getRelativeError() < 0.01) {
                compute_covariance_matrix(points, covarianceMatrix);
                choleskyDecomposition(covarianceMatrix);
            }

            if (slidingWindowStop.getRelativeError() < error)
                break;

//            std::cout << " $ " << interiorPoint.dot(objectiveFunction) << "\n";

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;

        return {minPoint, minPoint.dot(objectiveFunction)};
    }


    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    simulated_annealing_HMC(Spectrahedron &spectrahedron, Point &objectiveFunction,
                                             Parameters parameters, const NT error,
                                             const unsigned int maxSteps, Point &initial) {


        VT obj = objectiveFunction.getCoefficients();
        double normal = obj.norm();
        obj.normalize();
        Point objFunction = Point(obj);

        bool verbose = parameters.verbose;
        unsigned int walk_length = parameters.walk_steps;
        bool tillConvergence = maxSteps == 0;
        unsigned int step = 0;
        int dim = initial.dimension();
        std::list<Point> points;

        SlidingWindow slidingWindowStop(5 + sqrt(dim));

        Point interiorPoint = initial;

        double che_rad = maxDistanceFromPoint(spectrahedron, interiorPoint, parameters);

        double tempDescentFactor = 1 - 1 / (double) std::sqrt(dim);
        double temperature = che_rad;
        double temperature_threshold = 0.0001 / dim;

        double min = interiorPoint.dot(objFunction);
        Point minPoint = interiorPoint;
        Spectrahedron::BoundaryOracleBoltzmannHMCSettings<Point> settings;
        settings.first = true;
        settings.epsilon = 0.001;

        std::cout << " $ " << interiorPoint.dot(objectiveFunction) << "\n";

        // if fixed for debug Temp = che rad = 10 and initial point
//        che_rad = temperature = 10;
//        interiorPoint = initial; // for debug

        do {


            min_rand_point_generator_Boltzmann(spectrahedron, objFunction, interiorPoint, walk_length, che_rad, parameters,
                    temperature,  minPoint, min, settings);

            if (temperature > temperature_threshold) {
                temperature *= tempDescentFactor;
            }

            slidingWindowStop.push(interiorPoint.dot(objectiveFunction));


            if (slidingWindowStop.getRelativeError() < error)
                break;

//            std::cout << " $h " << interiorPoint.dot(objectiveFunction) << "\n";

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step << " steps" << std::endl;

        return {minPoint, minPoint.dot(objectiveFunction)};
    }
}
#endif //VOLESTI_SIMULATED_ANNEALING_SDP_H
