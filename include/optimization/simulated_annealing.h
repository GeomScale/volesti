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

    template<class Parameters, class Point, typename NT>
    std::pair<Point, NT>
    simulated_annealing(HPolytope<Point> &polytope, Point &objectiveFunction, Parameters parameters, const NT error,
                        const unsigned int maxSteps, Point &initial) {


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
        double min = interiorPoint.dot(objectiveFunction);
        Point minPoint = interiorPoint;

        std::cout << covarianceMatrix << "\n";
        do {
            temperature *= tempDescentFactor;

            min_rand_point_generator_Boltzmann(polytope, objectiveFunction, interiorPoint, walk_length, parameters, temperature, covarianceMatrix, minPoint, min);
            slidingWindow.push(min);

            std::cout << "========= step " << step  << " cost " << interiorPoint.dot(objectiveFunction) << " " << temperature << "=============\n";

            if (slidingWindow.getRelativeError() < error) {
                    Point temp = interiorPoint;
                    update_convariance_matrix(polytope, objectiveFunction, temperature, temp, parameters, rnum, walk_length, covarianceMatrix);
                    choleskyDecomposition(covarianceMatrix);
            }


            if (!polytope.is_in(interiorPoint)) {
                std::cout << "out\n";
                break;
            }

            step++;
        } while (step <= maxSteps || tillConvergence);

        STEPS = step;

        if (verbose) std::cout << "Ended at " << step<< " steps"  <<  std::endl;

        return {minPoint, min};
    }

}

#endif //VOLESTI_SIMULATED_ANNEALING_H
