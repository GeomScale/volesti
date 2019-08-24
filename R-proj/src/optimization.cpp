#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <chrono>
#include "truncated_exponential.h"
#include "lp_problem.h"
#include "sdp_problem.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "vars.h"

//' The function for linear proramming, i.e.
//'
//' min \eqn{cx}
//' s.t. \eqn{Ax\leq b}
//'
//' There are two available algorithms, the randomized cutting plane and the simulated annealing. The constraints of the linear program describe an H-polytope with \eqn{m} facets.
//'
//' @param P. An object from class Hpolytope.
//' @param objectiveFunction. A numerical vector with the coefficients of the objective function.
//' @param algorithm Optional. A string specifying which algorithm to use:
//' \itemize{
//' \item{\code{'RCP'} }{The randomized cutting plane method}
//' \item{\code{'RCPS'} }{The randomized cutting plane method with the implicit isotropization heuristic}
//' \item{\code{'SA'} }{The simulated annealing algorithm}
//' }
//' @param randomWalk Optional. A string specifying which random walk to use to use:
//' \itemize{
//' \item{\code{'BW'} }{The billiard walk, available only for \code{RCP}}
//' \item{\code{'CDHR'} }{Hit & Run with coordinate directions, available only for \code{RCP}}
//' \item{\code{'RDHR'} }{Hit & Run with random directions, available for \code{RCP}, \code{RCPS}, \code{SA}}
//' }
//' @param numMaxSteps Optional. The maximum number of iterations of the algorithm.
//' @param pointsNum Optional. The number of points to sample per iteration.
//' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{0.000001}.
//'
//'
//' @return A list with the minimum value and the values of x.
//' @examples
//' # calling RCP algorithm for a H-polytope (10d unit simplex)
//' P = volesti::GenSimplex(10, "H")
//' obj = c(1.80429, -8.46931, -1.68169, -1.71464, -1.95775, -4.24238, -7.19885, -1.64979, -5.9658, -1.189649)
//' volesti::lp_optimization(P, obj, algorithm = "SA", pointsNum=500)
//' @export
// [[Rcpp::export]]
Rcpp::List lp_optimization(Rcpp::Reference P,
                   Rcpp::NumericVector objectiveFunction,
                   Rcpp::Nullable<std::string> algorithm=R_NilValue,
                   Rcpp::Nullable<std::string> randomWalk=R_NilValue,
                   Rcpp::Nullable<unsigned int> numMaxSteps=R_NilValue,
                   Rcpp::Nullable<unsigned int> pointsNum=R_NilValue,
                   Rcpp::Nullable<double> error=R_NilValue)
{
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef optimization::lp_problem<Point, NT> lp_problem;
    typedef lp_problem::Algorithm Algorithm;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    /* CONSTANTS */
    //error in hit-and-run bisection of P
    const NT err = 0.0000000001;

    Algorithm algo;

    bool verbose = false,
            rand_only = false,
            round_only = false,
            file = false,
            round = false,
            NN = false,
            user_walk_len = false,
            linear_extensions = false,
            birk = false,
            rotate = false,
            ball_walk = false,
            ball_rad = false,
            experiments = true,
            annealing = false,
            Vpoly = false,
            Zono = false,
            cdhr = false,
            rdhr = true,
            exact_zono = false,
            gaussian_sam = false,
            billiard = false;

    lp_problem lp;
    int rnum, maxSteps;

    NT e = 1;
    NT distance = 0.0001;
    NT delta = -1.0;

    /** CHECK POLYTOPE **/

    Hpolytope HP;
    int dimension = Rcpp::as<Rcpp::Reference>(P).field("dimension");
    HP.init(dimension, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));

    /** BUILD LP **/

    VT obj(dimension);
    for (unsigned int j=0; j<dimension; j++){
        obj(j) = objectiveFunction[j];
    }

    lp = optimization::lp_problem<Point, NT>(HP, obj, optimization::Goal::minimize);

    /** CHECK ALGORITHM **/

    if(!algorithm.isNotNull() || Rcpp::as<std::string>(algorithm).compare(std::string("RCP"))==0){
        algo = Algorithm::RANDOMIZED_CUTTING_PLANE;
    }
    else if (Rcpp::as<std::string>(algorithm).compare(std::string("RCPS"))==0) {
        algo = Algorithm::RANDOMIZED_CUTTING_PLANE_SAMPLED_COVARIANCE_HEURISTIC;
    }
    else if (Rcpp::as<std::string>(algorithm).compare(std::string("SA"))==0) {
        algo = Algorithm::SIMULATED_ANNEALING_EFICIENT_COVARIANCE;
    }
    else {
        throw Rcpp::exception("Unknown method!");
    }


    /** CHECK RANDOM WALK **/

    if(!randomWalk.isNotNull() || Rcpp::as<std::string>(randomWalk).compare(std::string("CDHR"))==0){
        cdhr = true;
        rdhr = false;
    } else if (Rcpp::as<std::string>(randomWalk).compare(std::string("RDHR"))==0) {
        cdhr = false;
        rdhr = true;
        ball_walk = false;
    } else if (Rcpp::as<std::string>(randomWalk).compare(std::string("BW"))==0) {
        if (algo != Algorithm::RANDOMIZED_CUTTING_PLANE)
            throw Rcpp::exception("Can use billiard walk only with the randomized cutting plane method (RCP)!");

        algo = Algorithm::RANDOMIZED_CUTTING_PLANE_BILLIARD;
    } else {
        throw Rcpp::exception("Unknown walk type!");
    }

    /** CHECK MAX NUMBER OF STEPS **/

    if(!numMaxSteps.isNotNull()){
        maxSteps = 10*dimension;
    } else {
        maxSteps = Rcpp::as<NT>(numMaxSteps);
    }


    /** CHECK NUMBER OF POINTS **/

    if (pointsNum.isNotNull()) {
        rnum = Rcpp::as<unsigned int>(pointsNum);
    }
    else {
        switch (algo){
            case Algorithm::RANDOMIZED_CUTTING_PLANE:
            case Algorithm::RANDOMIZED_CUTTING_PLANE_SAMPLED_COVARIANCE_HEURISTIC:
            case Algorithm::SIMULATED_ANNEALING_EFICIENT_COVARIANCE:
                rnum = 1000 + dimension*dimension;
                break;
            case Algorithm::RANDOMIZED_CUTTING_PLANE_BILLIARD:
                rnum = 10;
                break;
        }
    }

    /** CHECK ERROR **/

    if(!error.isNotNull()){
        distance = 0.000001;
    } else {
        distance = Rcpp::as<NT>(error);
    }


    /**  SOLVE  **/

    /* RANDOM NUMBERS */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    vars<NT, RNGType> var(rnum, dimension, rnum, 1, err, e, 0, 0.0, 0, 0, rng,
                          urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);

    lp.solve(var, distance, maxSteps, algo);

    return Rcpp::List::create(Rcpp::Named("minvalue") = lp.solutionVal,
                                    Rcpp::Named("coordinates") = lp.solution.getCoefficients());
}


//' The function for semidefinite proramming, i.e.
//'
//' There are two available algorithms, the randomized cutting plane and the simulated annealing.
//'
//' @param S. An object from class Spectrahedron.
//' @param objectiveFunction. A numerical vector with the coefficients of the objective function.
//' @param algorithm Optional. A string specifying which algorithm to use:
//' \itemize{
//' \item{\code{'RCP'} }{The randomized cutting plane method}
//' \item{\code{'RCPS'} }{The randomized cutting plane method with the implicit isotropisation heuristic}
//' \item{\code{'SA'} }{The simulated annealing algorithm}
//' }
//' @param numMaxSteps Optional. The maximum number of iterations of the algorithm.
//' @param pointsNum Optional. The number of points to sample per iteration.
//' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{0.000001}.
//'
//'
//' @return A list with the minimum value and the values of x.
//' @export
// [[Rcpp::export]]
Rcpp::List sdp_optimization(Rcpp::Reference S,
                           Rcpp::NumericVector objectiveFunction,
                           Rcpp::Nullable<std::string> algorithm=R_NilValue,
                           Rcpp::Nullable<unsigned int> numMaxSteps=R_NilValue,
                           Rcpp::Nullable<unsigned int> pointsNum=R_NilValue,
                           Rcpp::Nullable<double> error=R_NilValue)
{
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef optimization::sdp_problem<Point> sdp_problem;
    typedef sdp_problem::Algorithm Algorithm;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    /* CONSTANTS */
    //error in hit-and-run bisection of P
    const NT err = 0.0000000001;

    Algorithm algo;

    bool verbose = false,
            rand_only = false,
            round_only = false,
            file = false,
            round = false,
            NN = false,
            user_walk_len = false,
            linear_extensions = false,
            birk = false,
            rotate = false,
            ball_walk = false,
            ball_rad = false,
            experiments = true,
            annealing = false,
            Vpoly = false,
            Zono = false,
            cdhr = false,
            rdhr = true,
            exact_zono = false,
            gaussian_sam = false,
            billiard = false;

    sdp_problem sdp;
    int rnum, maxSteps;

    NT e = 1;
    NT distance = 0.0001;
    NT delta = -1.0;

    /** CHECK SCPECTRAHEDRON **/

    Spectrahedron spectrahedron;
    std::vector<MT> matrices = Rcpp::as<Rcpp::Reference>(S).field("matrices");
    LMI lmi;
    lmi = LMI(matrices);
    spectrahedron = Spectrahedron(lmi);
    int dimension = Rcpp::as<Rcpp::Reference>(S).field("dimension");

    /** BUILD LP **/

    VT obj(dimension);
    for (unsigned int j=0; j<dimension; j++){
        obj(j) = objectiveFunction[j];
    }

    sdp = optimization::sdp_problem<Point>(spectrahedron, obj);

    /** CHECK ALGORITHM **/

    if(!algorithm.isNotNull() || Rcpp::as<std::string>(algorithm).compare(std::string("RCP"))==0){
        algo = Algorithm::RANDOMIZED_CUTTING_PLANE;
    }
    else if (Rcpp::as<std::string>(algorithm).compare(std::string("RCPS"))==0) {
        algo = Algorithm::RANDOMIZED_CUTTING_PLANE_COVARIANCE_MATRIX;
    }
    else if (Rcpp::as<std::string>(algorithm).compare(std::string("SA"))==0) {
        algo = Algorithm::SIMULATED_ANNEALING_EFICIENT_COVARIANCE;
    }
    else {
        throw Rcpp::exception("Unknown method!");
    }



    /** CHECK MAX NUMBER OF STEPS **/

    if(!numMaxSteps.isNotNull()){
        maxSteps = 10*dimension;
    } else {
        maxSteps = Rcpp::as<NT>(numMaxSteps);
    }


    /** CHECK NUMBER OF POINTS **/

    if (pointsNum.isNotNull()) {
        rnum = Rcpp::as<unsigned int>(pointsNum);
    }
    else {
        switch (algo){
            case Algorithm::RANDOMIZED_CUTTING_PLANE:
            case Algorithm::RANDOMIZED_CUTTING_PLANE_COVARIANCE_MATRIX:
            case Algorithm::SIMULATED_ANNEALING_EFICIENT_COVARIANCE:
                rnum = 1000 + dimension*dimension;
                break;
        }
    }

    /** CHECK ERROR **/

    if(!error.isNotNull()){
        distance = 0.000001;
    } else {
        distance = Rcpp::as<NT>(error);
    }


    /**  SOLVE  **/

    /* RANDOM NUMBERS */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    vars<NT, RNGType> var(rnum, dimension, rnum, 1, err, e, 0, 0.0, 0, 0, rng,
                          urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);

    sdp.solve(var, distance, maxSteps, algo);

    return Rcpp::List::create(Rcpp::Named("minvalue") = sdp.solution.second,
                              Rcpp::Named("coordinates") = sdp.solution.first.getCoefficients());

}