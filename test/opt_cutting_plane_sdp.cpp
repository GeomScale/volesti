
#include "Eigen/Eigen"

#define VOLESTI_DEBUG


#include <fstream>
#include "volume.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "solve_lp.h"

#include "cutting_plane_sdp.h"
#include "sdp_problem.h"
#include <chrono>
#include "Eigen"
#include <iomanip>

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef boost::mt19937 RNGType;
typedef optimization::sdp_problem<Point> sdp_problem;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

extern int STEPS; //TODO delete this

void printHelpMessage();

int main(const int argc, const char **argv) {

    // the object function is a vector

    //Deafault values
    int dimensinon, numOfExperinments = 1, walkLength = 10, numOfRandomPoints = 16, nsam = 100, numMaxSteps = 100;
    NT e = 1;
    bool useIsotropyMatrix = false;

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


    sdp_problem sdp;

    NT delta = -1.0, error = 0.2;
    NT distance = 0.0001;

    //parse command line input vars
    for (int i = 1; i < argc; ++i) {
        bool correct = false;


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
            sdp = sdp_problem(inp);
            inp.close();
            file = true;
            correct = true;
        }

        if (!strcmp(argv[i], "-lp")) {
            std::ifstream inp;
            inp.open(argv[++i], std::ios_base::in);
            sdp.transformFromLP(inp);
            inp.close();
//            sdp.print();
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
            useIsotropyMatrix = true;
            correct = true;
        }

        if (!strcmp(argv[i], "-r")) {
            numOfRandomPoints = atoi(argv[++i]);
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

    // Setup the parameters
    vars<NT, RNGType> var(numOfRandomPoints, dimensinon, walkLength, 1, err, e, 0, 0.0, 0, 0, rng,
                          urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);

    //RUN EXPERIMENTS
    std::vector<NT> results;
    std::vector<double> times;
    std::vector<int> steps;

    for (unsigned int i = 0; i < numOfExperinments; ++i) {
//        std::cout << "Experiment " << i + 1 << std::endl;
//        auto t1 = std::chrono::steady_clock::now();
//
//        lp.solve(useIsotropyMatrix, var, distance, numMaxSteps);
//
//        auto t2 = std::chrono::steady_clock::now();
//
//        std::cout << std::fixed;
//        lp.printSolution();
//
//        if ( std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() < 10000 ) {
//            std::cout << "Computed at " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msecs" << std::endl << std::endl;
//        }
//        else {
//            std::cout << "Computed at " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " secs" << std::endl << std::endl;
//        }
//
//        results.push_back(lp.solutionVal);
//        steps.push_back(STEPS);
//        times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
        sdp.solve(var, distance, numMaxSteps);
        std::cout << std::fixed;
        std::cout << std::setprecision(3);
        sdp.printSolution();
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

//    auto t1 = std::chrono::steady_clock::now();
//    double exact = solveWithLPSolve(lp);
//    auto t2 = std::chrono::steady_clock::now();
//    std::cout << "\nStatistics\n" <<
//              "Average result: " << average << "\n"<<
//              "Average time: " << avg_time << "\n" <<
//              "Average # Steps: " << steps_avg << "\n" <<
//              "Variance: " << variance << "\n" <<
//              "Standard deviation: " << std_dev << "\n" <<
//              "Coefficient of variation: " << abs(std_dev / average)  << "\n" <<
//              "Exact Solution: " << exact  << "\n" <<
//              "lpsolve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " \n"
//                                                                                                            "Relative error: " << abs((exact - average) / exact) << "\n";
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
