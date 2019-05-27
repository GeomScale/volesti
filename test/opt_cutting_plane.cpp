
#include "Eigen/Eigen"

#define VOLESTI_DEBUG


#include <fstream>
#include "volume.h"
#include "sample_only.h"
#include "exact_vols.h"
#include "cutting_plane.h"
#include "solve_lp.h"
#include <chrono>
#include "Eigen"

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef boost::mt19937 RNGType;
typedef HPolytope<Point> Hpolytope;
typedef std::pair<Point, NT> Result;

typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

void printHelpMessage();

bool
readFromFile(const char *const *argv, bool verbose, HPolytope<point<Cartesian<double>::Self>> &HP, int &n, bool &file,
             int &i, VT& objectFunction);


bool loadProgramFromStream(std::istream &is, HPolytope<Point> &HP, VT& objectFunction);
NT solveWithLPSolve(HPolytope<Point>& HP, VT objectFunction);

int main(const int argc, const char **argv) {

    // the object function is a vector

    VT objectFunction;
    //Deafault values
    int dimensinon, numOfExperinments = 1, walkLength = 10, numOfRandomPoints = 16, nsam = 100, numMaxSteps = 100;
    NT e = 1;
    bool uselpSolve = false;
    bool useIsotropyMatrix = true;

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


    Hpolytope HP;

    NT delta = -1.0, error = 0.2;
    NT distance = 0.0001;

    //parse command line input vars
    for (int i = 1; i < argc; ++i) {
        bool correct = false;

        if (!strcmp(argv[i], "-lpsolve")) {
            uselpSolve = true;
            correct = true;
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
            readFromFile(argv, verbose, HP, dimensinon, file, i, objectFunction);
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
        if (!strcmp(argv[i], "-noISO")) {
            useIsotropyMatrix = false;
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



    if (uselpSolve) {
        std::cout << "Using lp_solve" << std::endl;

        auto t1 = std::chrono::steady_clock::now();


        NT min = solveWithLPSolve(HP, objectFunction);

        auto t2 = std::chrono::steady_clock::now();

        std::cout << "Computed at: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msecs" << std::endl;

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


    //RUN EXPERIMENTS
    std::vector<double> times;
    std::vector<Result> results;
    NT sum = NT(0);
    std::cout.precision(7);


    if (verbose && HP.num_of_hyperplanes() < 100) {
        std::cout << "Input polytope is of dimension: " << dimensinon << std::endl;
        HP.print();
    }

    for (unsigned int i = 0; i < numOfExperinments; ++i) {
        std::cout << "Experiment " << i + 1 << std::endl;

        auto t1 = std::chrono::steady_clock::now();

        // Setup the parameters
        vars<NT, RNGType> var(numOfRandomPoints, dimensinon, walkLength, 1, err, e, 0, 0.0, 0, 0, rng,
                              urdist, urdist1, delta, verbose, rand_only, round, NN, birk, ball_walk, cdhr, rdhr);



        Result result;

        if (useIsotropyMatrix)
            result = optimization::cutting_plane_method_isotropic(objectFunction, HP, var, distance, numMaxSteps);
        else
            result = optimization::cutting_plane_method(objectFunction, HP, var, distance, numMaxSteps);

        auto t2 = std::chrono::steady_clock::now();


        if ( std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() < 10000 ) {
            std::cout << "Min value is: " << result.second << std::endl <<
                      "coords: ";
            result.first.print();
            std::cout << "Computed at " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " msecs" << std::endl << std::endl;
        }
        else {
            std::cout << "Min value is: " << result.second << std::endl <<
                      "coords: ";
            result.first.print();
            std::cout << "Computed at " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " secs" << std::endl << std::endl;
        }

//        results.push_back(result);
//        sum += result.second;
    }


    return 0;
}



bool
readFromFile(const char *const *argv, bool verbose, HPolytope<point<Cartesian<double>::Self>> &HP, int &n, bool &file,
             int &i, VT& objectFunction) {
    file = true;
    std::cout << "Reading input from file..." << std::endl;
    std::ifstream inp;
    inp.open(argv[++i], std::ios_base::in);
    bool retval = loadProgramFromStream(inp, HP, objectFunction);
    inp.close();
    n = HP.dimension();
//    HP.print();
    return retval;
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

bool loadProgramFromStream(std::istream &is, HPolytope<Point> &HP, VT& objectFunction){

    std::string line;
    std::string::size_type sz;

    //read dimension
    if (std::getline(is, line, '\n').eof())
        return false;

    int dim = std::stoi(line);


    //read object function
    if (std::getline(is, line, '\n').eof())
        return false;

    objectFunction.resize(dim);

    NT num = std::stod(line, &sz);
    objectFunction(0) = num;

    for (int j=2 ; j<=dim  ; j++) {
        line = line.substr(sz);
        num = std::stod(line, &sz);
        objectFunction(j-1) = num;
    }


    if (std::getline(is, line, '\n').eof())
        return false;

    //read number of constraints
    int constraintsNum = std::stoi(line, &sz);

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    VT b;
    MT A;

    A.resize(constraintsNum, dim);
    b.resize(constraintsNum);

    // for each constraint
    for (int i=1 ; i<=constraintsNum ; i++) {

        if (std::getline(is, line, '\n').eof())
            return false;

        // read first line of A

        NT num = std::stod(line, &sz);
        A(i - 1, 0) = num;

        for (int j=2 ; j<=dim  ; j++) {
            line = line.substr(sz);
            num = std::stod(line, &sz);
            A(i - 1, j - 1) = num;
        }

        //read first row of b
        line = line.substr(sz);
        num = std::stod(line, &sz);
        b(i - 1) = num;
    }

    HP.init(dim, A, b);
    return true;
}

NT solveWithLPSolve(HPolytope<Point>& HP, VT objectFunction) {
    lprec *lp;
    unsigned int dim = HP.dimension();

    REAL row[1 + dim]; /* must be 1 more then number of columns ! */

    /* Create a new LP model */
    lp = make_lp(0, dim);

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    VT b = HP.get_vec();
    MT A = HP.get_mat();


    for (int j=1 ; j<=dim ; j++)
        row[j] = objectFunction(j-1); //j must start at 1

    set_obj_fn(lp, row);

    set_add_rowmode(lp, TRUE);

    for (int i=0 ; i<A.rows() ; i++) {
        for (int j=1 ; j<=dim ; j++)
            row[j] = A(i, j-1); //j must start at 1

        add_constraint(lp, row, LE, b(i)); /* constructs the row: +v_1 +2 v_2 >= 3 */
    }

    set_add_rowmode(lp, FALSE);
    set_minim(lp);

    solve(lp);
    NT ret = get_objective(lp);
    delete_lp(lp);
    return ret;
}