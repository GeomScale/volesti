
#include "Eigen/Eigen"

#define VOLESTI_DEBUG


#include <fstream>
#include "volume.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "solve_lp.h"
#include <chrono>
#include "Eigen"
#include "lp_problem.h"
#include "lp_generator.h"
#include "interior_point.h"
#include <iomanip>

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef boost::mt19937 RNGType;
typedef HPolytope<Point> Hpolytope;
typedef optimization::lp_problem<Point, NT> lp_problem;
typedef lp_problem::Algorithm Algorithm;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

extern int STEPS; //TODO delete this

void printHelpMessage();
NT solveWithLPSolve(lp_problem);

int main(const int argc, const char **argv) {

    // the object function is a vector

    //Deafault values
    int dimensinon, numOfExperinments = 1, walkLength = 10, numOfRandomPoints = 16, nsam = 100, numMaxSteps = 100;
    NT e = 1;
    bool uselpSolve = false;
    Algorithm algorithm = Algorithm::RANDOMIZED_CUTTING_PLANE;

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
            rdhr = true, // for hit and run
            exact_zono = false,
            gaussian_sam = false;


    lp_problem lp;

    NT delta = -1.0, error = 0.2;
    NT distance = 0.0001;

    //parse command line input vars
    for (int i = 1; i < argc; ++i) {
        bool correct = false;

        if (!strcmp(argv[i], "-lpsolve")) {
            uselpSolve = true;
            correct = true;
        }

        if (!strcmp(argv[i], "-cdhr")) {
            cdhr = true;
            correct = true;
            rdhr = false; // for hit and run
        }
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printHelpMessage();
            return 0;
        }

        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
            std::cout << "Verbose mode\n";
            correct = true;
        }

        if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--file")) {
            std::ifstream inp;
            inp.open(argv[++i], std::ios_base::in);
            lp = optimization::lp_problem<Point, NT>(inp);
            inp.close();
            file = true;
            correct = true;
        }

        if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--error")) {
            e = atof(argv[++i]);
            distance = e;
            correct = true;
        }
        if (!strcmp(argv[i], "-w") || !strcmp(argv[i], "--walk_len")) {
            walkLength = atoi(argv[++i]);
            user_walk_len = true;
            correct = true;
        }
        if (!strcmp(argv[i], "-ISO")) {
            algorithm = Algorithm::RANDOMIZED_CUTTING_PLANE_SAMPLED_COVARIANCE_HEURISTIC;
            correct = true;
        }

        if (!strcmp(argv[i], "-r")) {
            numOfRandomPoints = atoi(argv[++i]);
            correct = true;
        }

        if (!strcmp(argv[i], "-deterministic")) {
            algorithm = Algorithm::DETERMINISTIC_CUTTING_PLANE_CHEBYSHEV_CENTER;
            correct = true;
        }

        if (!strcmp(argv[i], "-simulatedannealing")) {
            algorithm = Algorithm::SIMULATED_ANNEALING;
            correct = true;
        }

        if (!strcmp(argv[i], "-exp")) {
            numOfExperinments = atoi(argv[++i]);
            correct = true;
        }

        if (!strcmp(argv[i], "-k")) {
            numMaxSteps = atoi(argv[++i]);
            correct = true;
        }

        if (!correct) {
            std::cerr << "unknown parameters \'" << argv[i] <<
                      "\', try " << argv[0] << " --help" << std::endl;
            exit(-2);
        }

    }


    if (uselpSolve) {
        std::cout << "Using lp_solve" << std::endl;

        auto t1 = std::chrono::steady_clock::now();


        NT min = solveWithLPSolve(lp);

        auto t2 = std::chrono::steady_clock::now();

        std::cout << "Min is " << min << std::endl << "Computed at: " <<  std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msecs" << std::endl;
        return 0;
    }

    /* CONSTANTS */
    //error in hit-and-run bisection of P
    const NT err = 0.0000000001;


    // If no file specified
    if (!file) {
        std::cout << "You must specify a file - type -h for help" << std::endl;
        exit(-2);
    }


    /* RANDOM NUMBERS */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::normal_distribution<> rdist(0, 1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

//    for (int j = 0; j < 10; ++j) {
//        std::cout << texp(e, numOfRandomPoints, walkLength, rng) << "\n";
//    }
//
//    return 0;
    // Setup the parameters
    vars<NT, RNGType> var(numOfRandomPoints, dimensinon, walkLength, 1, err, e, 0, 0.0, 0, 0, rng,
                          urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);

    //RUN EXPERIMENTS
    std::vector<NT> results;
    std::vector<double> times;
    std::vector<int> steps;

    for (unsigned int i = 0; i < numOfExperinments; ++i) {
        std::cout << "Experiment " << i + 1 << std::endl;
        auto t1 = std::chrono::steady_clock::now();

        lp.solve(var, distance, numMaxSteps, algorithm);

        auto t2 = std::chrono::steady_clock::now();

        std::cout << std::fixed;
        lp.printSolution();

        if ( std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() < 10000 ) {
            std::cout << "Computed at " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msecs" << std::endl << std::endl;
        }
        else {
            std::cout << "Computed at " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " secs" << std::endl << std::endl;
        }

        results.push_back(lp.solutionVal);
        steps.push_back(STEPS);
        times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
    }

    double sum = 0;

    for (auto x : results)
        sum += x;

    double average = sum / (double) results.size();
    double std_dev=0;
    double variance = 0;

    for(auto x : results){
        variance += std::pow(x - average,2);
    }

    variance /= (double)results.size();
    std_dev = std::sqrt(variance);

    double avg_time = 0;

    for (auto x : times)
        avg_time += x;

    avg_time /= (double) times.size();

    int steps_avg = 0;

    for (auto x : steps)
        steps_avg += x;

    steps_avg = steps_avg / steps.size();

    std::cout << std::fixed;
    std::cout << std::setprecision(3);

    auto t1 = std::chrono::steady_clock::now();
    double exact = solveWithLPSolve(lp);
    auto t2 = std::chrono::steady_clock::now();
    std::cout << "\nStatistics\n" <<
        "Average result: " << average << "\n"<<
         "Average time: " << avg_time << "\n" <<
         "Average # Steps: " << steps_avg << "\n" <<
        "Variance: " << variance << "\n" <<
        "Standard deviation: " << std_dev << "\n" <<
        "Coefficient of variation: " << abs(std_dev / average)  << "\n" <<
        "Exact Solution: " << exact  << "\n" <<
        "lpsolve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " \n"
        "Relative error: " << abs((exact - average) / exact) << "\n";
    return 0;
}


void printHelpMessage() {
    std::cerr <<
              "Usage: The constraints are passed in a file, as in vol.cpp while the object function is declared in main\n" <<
              "-v, --verbose \n" <<
              "-r [#num]: the number of points to sample at each k (default 16)\n" <<
              "-k [#num]: the number of maximum iterations (default 100) - if 0 runs until convergence\n" <<
              "-f, --file [filename] The file must be of the following format:\n" <<
              "\t #dimension\n" <<
              "\t object function\n" <<
              "\t #num_of_constaints\n" <<
              "\t constraints\n" <<
              "-e, --error epsilon : the goal error of approximation\n" <<
              "-w, --walk_len [walk_len] : the random walk length (default 10)\n" <<
              "-exp [#exps] : number of experiments (default 1)\n" <<
              "-d [#distance] : stop if successive estimations are less than (default 0.000001)\n" <<
              "-noISO: don't use isotropy matrix" <<
              std::endl;
}

NT solveWithLPSolve(lp_problem lp) {
    lprec *_lp;
    unsigned int dim = lp.polytope.dimension();

    REAL row[1 + dim]; /* must be 1 more then number of columns ! */

    /* Create a new LP model */
    _lp = make_lp(0, dim);

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    VT b = lp.polytope.get_vec();
    MT A = lp.polytope.get_mat();


    for (int j=1 ; j<=dim ; j++)
        row[j] = lp.objectiveFunction(j-1); //j must start at 1

    set_obj_fn(_lp, row);
    set_verbose(_lp, 2);
    set_add_rowmode(_lp, TRUE);

    for (int i=0 ; i<A.rows() ; i++) {
        for (int j=1 ; j<=dim ; j++)
            row[j] = A(i, j-1); //j must start at 1

        add_constraint(_lp, row, LE, b(i)); /* constructs the row: +v_1 +2 v_2 >= 3 */
    }

    set_add_rowmode(_lp, FALSE);
    set_minim(_lp);

    for (int j=1 ; j<=dim ; j++)
        set_bounds(_lp, j, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max());

    solve(_lp);
    NT ret = get_objective(_lp);

    REAL solution[dim];
    get_variables(_lp, solution);

//    std::cout << "Optimal solution: " << std::endl;
//    for (int i=0; i < dim; i++)
//        std::cout << "x[" << i+1 << "] = " << solution[i] << std::endl;
//        std::cout << " , " <<  solution[i];

    delete_lp(_lp);
    return ret;
}