//
// Created by panagiotis on 31/7/2019.
//

#ifndef VOLESTI_SIMULATED_ANNEALING_H
#define VOLESTI_SIMULATED_ANNEALING_H

#include "polytopes.h"
#include "samplers.h"
#include "heuristics.h"
#include "Eigen"

namespace optimization {

    void choleskyDecomposition(MT& covarianceMatrix) {
        Eigen::LLT<MT> lltOfA(covarianceMatrix); // compute the Cholesky decomposition of A
        covarianceMatrix = lltOfA.matrixL();
    }

    template <class Point>
    void compute_covariance_matrix(std::list<Point> points, MT& covarianceMatrix) {
        int dim = points.front().dimension();
        covarianceMatrix.setZero(dim, dim);

        double factor = 1 / (double) points.size();
        Point sum(dim);

        for (auto p : points) {
            sum = sum + p;
        }
        sum = sum * factor;

        VT sumCoeffs = sum.getCoefficients();

        for (auto p : points) {
            VT coeffs = p.getCoefficients();
            covarianceMatrix += (coeffs * coeffs.transpose());
        }

        covarianceMatrix *= factor;
        covarianceMatrix -= sumCoeffs * sumCoeffs.transpose();
    }

    template<class Polytope, class Parameters, class Point>
    void init_convariance_matrix_uniform(Polytope &P, const Point& p, Parameters& vars, unsigned int r_num, unsigned int walk_len, MT& covariance_matrix) {
        std::list<Point> points;
        Point temp = p;

        for (int i=0 ; i<r_num ; i++) {
            for (int k=0 ; k<walk_len ; k++) {
                hit_and_run(temp, P, vars);
            }

            points.push_back(temp);
        }

        compute_covariance_matrix(points, covariance_matrix);
    }


    template<class Polytope, class Parameters, class Point>
    void update_convariance_matrix(Polytope &P, Point& c, double temperature, Point& p, Parameters& vars, unsigned int r_num, unsigned int walk_len, MT& covariance_matrix) {
        std::list<Point> points;
        rand_point_generator_Boltzmann(P, c, p, r_num, walk_len, vars, temperature, covariance_matrix, points);

        compute_covariance_matrix(points, covariance_matrix);
    }

    template <class Polytope, class Point, class Parameters>
    void min_hit_and_run_Boltzmann(Point &p,
                                   Polytope &P,
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

        std::pair <NT, NT> dbpair = P.line_intersect(p, l);
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


    template <class Polytope, class Point, class Parameters>
    void min_hit_and_run_Boltzmann(Point &p,
                               Polytope &P,
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

        std::pair <NT, NT> dbpair = P.line_intersect(p, l);
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

    template <class Polytope, class Point, class Parameters>
    void min_hit_and_run_Boltzmann(Point &p,
                                   Polytope &P,
                                   Parameters &var,
                                   Point& BoltzmannDirection,
                                   double BoltzmannParameter,
                                   const MT& choleskyDecomp,
                                   Point &minPoint,
                                   double& minValue) {
        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;
        unsigned int n =p.dimension();
        RNGType &rng = var.rng;

        Point l = get_direction<RNGType, Point, NT>(n, choleskyDecomp);

        std::pair <NT, NT> dbpair = P.line_intersect(p, l);
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

    template <class Polytope, class Point, class Parameters>
    void min_hit_and_run_Boltzmann_set_direction(Point &p,
                                   Polytope &P,
                                   Parameters &var,
                                   Point& BoltzmannDirection,
                                   double BoltzmannParameter,
                                   Point& direction,
                                   Point &minPoint,
                                   double& minValue) {
        typedef typename Parameters::RNGType RNGType;
        typedef typename Point::FT NT;
        RNGType &rng = var.rng;
        std::pair <NT, NT> dbpair = P.line_intersect(p, direction);
        NT min_plus = dbpair.first;
        NT max_minus = dbpair.second;
        Point b1 = (min_plus * direction) + p;
        Point b2 = (max_minus * direction) + p;
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

        p = (lambda * direction) + p;
    }

    template<class Polytope, class Parameters, class Point>
    void min_rand_point_generator_Boltzmann(Polytope &P,
                                        Point& c,
                                        Point &p,   // a point to start
                                        unsigned int walk_length,
                                        Parameters &var,
                                        double temperature,
                                        const MT& covariance_matrix,
                                        Point &minPoint,
                                        double& minValue) {

        // begin sampling
        for (unsigned int i = 1; i <= walk_length ; ++i) {
            min_hit_and_run_Boltzmann(p, P, var, c, temperature, covariance_matrix, minPoint, minValue);
        }
    }

    template<class Polytope, class Parameters, class Point>
    void min_rand_point_generator_Boltzmann(Polytope &P,
                                            Point& c,
                                            Point &p,   // a point to start
                                            unsigned int walk_length,
                                            Parameters &var,
                                            double temperature,
                                            const MT& covariance_matrix,
                                            Point &minPoint,
                                            double& minValue,
                                            std::list<Point>& points,
                                            double& avgMinPerPhase) {

        // begin sampling
        for (unsigned int i = 0; i < walk_length ; ++i) {
            min_hit_and_run_Boltzmann(p, P, var, c, temperature, covariance_matrix, minPoint, minValue, avgMinPerPhase);
            points.push_back(p);
        }
    }

    template<class Polytope, class Parameters, class Point>
    void min_rand_point_generator_Boltzmann_window(Polytope &P,
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

        min_hit_and_run_Boltzmann(p, P, var, c, temperature, covariance_matrix, minPoint, minValue, avgMinPerPhase);
        slidingWindow.push(avgMinPerPhase);
        points.push_back(p);

        while (slidingWindow.getRelativeError() > 0.00001) {
            min_hit_and_run_Boltzmann(p, P, var, c, temperature, covariance_matrix, minPoint, minValue, avgMinPerPhase);
            count++;
            slidingWindow.push(avgMinPerPhase / (double) count);
            if (count % pickEverySteps == 0) points.push_back(p);
        }

//        std::cout << "-"<<count ;
    }


    template<class Polytope, class Parameters, class Point>
    void min_rand_point_generator_Boltzmann_set_directions(Polytope &P,
                                            Point& c,
                                            Point &p,   // a point to start
                                            unsigned int walk_length,
                                            Parameters &var,
                                            double temperature,
                                            Point &minPoint,
                                            double& minValue,
                                            std::vector<Point>& previous_points,
                                            std::vector<Point>& points) {
        typedef typename Parameters::RNGType RNGType;
        RNGType &rng = var.rng;
        int pointsNum = previous_points.size();
        boost::random::uniform_int_distribution<> uniformIntDistribution(0,pointsNum-1);
        Point previousP = p;

        // begin sampling
        // walk_length = pointsNum
        for (unsigned int i = 0; i < walk_length ; ++i) {
            int index = uniformIntDistribution(rng);
            Point direction = previous_points[index] - previousP;
            min_hit_and_run_Boltzmann_set_direction(p, P, var, c, temperature, direction, minPoint, minValue);
            points[i] = p;
        }
    }

    template<class Polytope, class Parameters, class Point>
    void min_rand_point_generator_Boltzmann(Polytope &P,
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
            min_hit_and_run_Boltzmann(p, P, var, c, temperature, minPoint, minValue);
            points.push_back(p);
        }
    }

    template<class Polytope, class Parameters, class Point>
    void min_rand_point_generator_Boltzmann(Polytope &P,
                                            Point& c,
                                            Point &p,   // a point to start
                                            unsigned int walk_length,
                                            Parameters &var,
                                            double temperature,
                                            Point &minPoint,
                                            double& minValue,
                                            std::vector<Point>& points) {

        // begin sampling
        for (unsigned int i = 0 ; i < walk_length ; ++i) {
            min_hit_and_run_Boltzmann(p, P, var, c, temperature, minPoint, minValue);
            points[i] = p;
        }
    }

    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    simulated_annealing(HPolytope<Point> &polytope, Point &objectiveFunction, Parameters parameters, const NT error,
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

        SlidingWindow slidingWindow(3);
        Point interiorPoint = initial;
        std::pair<double, double> radii = polytope.farthest_closest_distances(initial);
        MT covarianceMatrix;
        double temperature;
        bool updatedCovarianceMatrix = false;

        init_convariance_matrix_uniform(polytope, interiorPoint, parameters, rnum, walk_length, covarianceMatrix);
        choleskyDecomposition(covarianceMatrix);

        double tempDescentFactor = 1 - 1/(double) std::sqrt(dim);
        temperature = radii.first;
        double min = interiorPoint.dot(objFunction);
        Point minPoint = interiorPoint;

        do {
            temperature *= tempDescentFactor;

            min_rand_point_generator_Boltzmann(polytope, objFunction, interiorPoint, walk_length, parameters, temperature, covarianceMatrix, minPoint, min);
            slidingWindow.push(min);

//            std::cout << "========= step " << step  << " cost " << interiorPoint.dot(objFunction) << " " << temperature << "=============\n";
//            std::cout << interiorPoint.dot(objectiveFunction) << ", ";

            if (slidingWindow.getRelativeError() < error) {
                    Point temp = interiorPoint;
                    update_convariance_matrix(polytope, objFunction, temperature, temp, parameters, rnum, walk_length, covarianceMatrix);
                    choleskyDecomposition(covarianceMatrix);
            }


//            if (!polytope.is_in(interiorPoint)) {
//                std::cout << "out\n";
//                break;
//            }

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step<< " steps"  <<  std::endl;

        return {minPoint, minPoint.dot(objectiveFunction)};
    }

    template<class Point>
    void getArithmeticMean(std::vector<Point> points, Point& mean) {
        int dim = points[0].dimension();
        mean = Point(dim);

        for (auto p : points)
            mean = mean + p;

        mean = mean / points.size();
    }

    template<class Point>
    void getArithmeticMean(std::list<Point> points, Point& mean) {
        int dim = points.front().dimension();
        mean = Point(dim);

        for (auto p : points)
            mean = mean + p;

        mean = mean / points.size();
    }

    //the covariance matrix is computed by the intermediate points of the random walk
    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    simulated_annealing_efficient_covariance(HPolytope<Point> &polytope, Point &objectiveFunction, Parameters parameters, const NT error,
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
        std::pair<double, double> radii = polytope.farthest_closest_distances(initial);
        MT covarianceMatrix;

        double tempDescentFactor = 1 - 1/(double) std::sqrt(dim);
        double temperature = radii.first;
        double temperature_threshold = 0.000001 / dim;
        double avgMinPerPhase;

        double min = interiorPoint.dot(objFunction);
        Point minPoint = interiorPoint;

        min_rand_point_generator_Boltzmann(polytope, objFunction, interiorPoint, walk_length, parameters, temperature, minPoint, min, points);

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
            min_rand_point_generator_Boltzmann_window(polytope, objFunction, interiorPoint, error, parameters, temperature, covarianceMatrix, minPoint, min, points, avgMinPerPhase);

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

        if (verbose) std::cout << "Ended at " << step<< " steps"  <<  std::endl;

        return {minPoint, minPoint.dot(objectiveFunction)};
    }



    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    simulated_annealing_set_directions(HPolytope<Point> &polytope, Point &objectiveFunction, Parameters parameters, const NT error,
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
        std::vector<Point> points(walk_length);
        std::vector<Point> previous_points(walk_length);

        SlidingWindow slidingWindow(3);
        Point interiorPoint = initial;
        std::pair<double, double> radii = polytope.farthest_closest_distances(initial);
//        MT covarianceMatrix;
        double temperature;
//        bool updatedCovarianceMatrix = false;

//        init_convariance_matrix_uniform(polytope, interiorPoint, parameters, rnum, walk_length, covarianceMatrix);
//        choleskyDecomposition(covarianceMatrix);

        double tempDescentFactor = 1 - 1/(double) std::sqrt(dim);
        temperature = radii.first; //0.01
        double min = interiorPoint.dot(objFunction);
        Point minPoint = interiorPoint;

        min_rand_point_generator_Boltzmann(polytope, objFunction, interiorPoint, walk_length, parameters, temperature, minPoint, min, points);
//        compute_covariance_matrix(points, covarianceMatrix);
//        choleskyDecomposition(covarianceMatrix);

        double temperature_threshold = error / dim;

        do {
            if (temperature > temperature_threshold)
                temperature *= tempDescentFactor;

            previous_points = points;
            getArithmeticMean(previous_points, interiorPoint);
            min_rand_point_generator_Boltzmann_set_directions(polytope, objFunction, interiorPoint, walk_length, parameters, temperature, minPoint, min, previous_points, points);
            slidingWindow.push(min);

//            std::cout << "========= step " << step  << " cost " << interiorPoint.dot(objFunction) << " " << temperature << "=============\n";

//            if (slidingWindow.getRelativeError() < 0.01) {
//                compute_covariance_matrix(points, covarianceMatrix);
//                choleskyDecomposition(covarianceMatrix);
//            }

//            if (!polytope.is_in(interiorPoint)) {
//                std::cout << "out\n";
//                break;
//            }
//            std::cout << interiorPoint.dot(objectiveFunction) << ", ";

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step<< " steps"  <<  std::endl;

        return {minPoint, minPoint.dot(objectiveFunction)};
    }
}

#endif //VOLESTI_SIMULATED_ANNEALING_H
