// Copyright 2012-2013 National and Kapodistrian University of Athens, Greece.
//
// This file is part of RandGeom.
//
// RandGeom is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// RandGeom is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// See the file COPYING.LESSER for the text of the GNU Lesser General
// Public License.  If you did not receive this file along with RandGeom,
// see <http://www.gnu.org/licenses/>.
//
// Developer: Vissarion Fisikopoulos

#include <vol_rand.h>
#include <random>
#include <string>
#include <list>
#include <chrono>
//#include <proc/readproc.h>
class Timer {
public:
	Timer() { start_time = std::chrono::high_resolution_clock::now(); }

	void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}

	double end() {
		return elapsed_seconds();
	}

	double elapsed_seconds() {
		auto end_time = std::chrono::high_resolution_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
		return elapsed.count();
	}
private:
	decltype(std::chrono::high_resolution_clock::now()) start_time;
};

//////////////////////////////////////////////////////////
/**** MAIN *****/
//////////////////////////////////////////////////////////

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// Approximating the volume of a convex polytope or body
// can also be used for integration of concave functions.
// The user should provide the appropriate membership
// oracles.

int boundary_main(stdHPolytope<double>& P, double epsilon, int num_query_points, vars& var, bool use_jl) {
	Point chebPoint = P.create_point_representation();
	if (use_jl) {
		P.create_ann_jl_ds();
	}
	else {
		P.create_ann_ds();
	}

	std::list<Point> randPoints; //ds for storing rand points

	std::cout <<"Creating points" << std::endl;
	rand_point_generator(P, chebPoint, 2*num_query_points, var.walk_steps, randPoints, var);
	std::cout <<"Created points" << std::endl;

	auto it = randPoints.begin();
	double line_intersect_time = 0;
	double boundary_time = 0;
	std::vector<int> avgSteps;
	bool succeeded = false;
	double maxDist = -1;
	double minDist = 100000000;
	double avgDist = 0;
	int failed =0 ;
	for (int i=0; i<num_query_points; i++) {
		std::cout << ((double)i/num_query_points*100)<<"% completed " << std::endl;
		Point p = (*it);
		++it;
		Vector v = (*it) - CGAL::ORIGIN;
		//std::vector<double> v_v;
		//for (int i=0; i<P.dimension()-1; i++) {
		//	v_v.push_back(0);
		//}
		//v_v.push_back(2);
		//Vector v(P.dimension(), v_v.begin(), v_v.end());
		++it;
		Ray ray(p, v);
		//std::cout <<  "Ray source " << p << "\nRay direction " << v << std::endl;
		Timer timer1;
		auto range = P.line_intersect(p, v);
		double line_intersect_time_tmp = timer1.elapsed_seconds();
		line_intersect_time += line_intersect_time_tmp;

		Timer timer2;
		int numberOfSteps = 0;
		Point p2 = P.compute_boundary_intersection(ray, &numberOfSteps, &succeeded, epsilon, use_jl);
		double oracle_time = timer2.elapsed_seconds();
		boundary_time += oracle_time;
		if (succeeded) {
			avgSteps.push_back(numberOfSteps);
			double dist = std::sqrt(((p2-CGAL::ORIGIN) - (range.first-CGAL::ORIGIN)).squared_length()) ;
			if (dist<minDist) {
				minDist = dist;
			}
			if (dist>maxDist) {
				maxDist = dist;
				std::cout << "Changed max dist to " << maxDist << std::endl;
			}

			avgDist += dist;
			if (epsilon<=std::sqrt(((p2-CGAL::ORIGIN) - (range.first-CGAL::ORIGIN)).squared_length()) ) {

				//std::cout << "Distance to first " << std::sqrt(((p2-CGAL::ORIGIN) - (range.first-CGAL::ORIGIN)).squared_length()) << std::endl;
				//std::cout << "bp: " << p2 << std::endl;
				//std::cout << "ep: " << range.first << std::endl;
			}
		}
		else {
			failed++;
		}
	}

	std::cout << "Time line_intersect (" << line_intersect_time << "s) : " << std::endl;
	std::cout << "Time boundary oracle (" << boundary_time << "s) :" << std::endl;
	double sum = 0;
	for (int i=0; i<avgSteps.size(); i++) {
		sum += avgSteps[i];
	}
	std::cout << "Avg steps " << sum / avgSteps.size() << "\n";
	std::cout << "Failed attempts " << failed << std::endl;
	std::cout << "Min dist: " << minDist << std::endl;
	std::cout << "Max dist: " << maxDist << std::endl;
	std::cout << "Avg dist:  " << (double)avgDist/(num_query_points-failed) << std::endl;
}

int membership_main(stdHPolytope<double>& P, double epsilon, int num_query_points, vars& var, bool use_jl=true) {
	Point chebPoint = P.create_point_representation();

	Timer ann_create;
	if (use_jl) {
		P.create_ann_jl_ds();
	}
	else {
		P.create_ann_ds();
	}
	double ann_create_time = ann_create.elapsed_seconds();
	std::cout << "ANN " << (use_jl?"(with jl) ": " ") << "took " << ann_create_time << "s to build" << std::endl;

	std::list<Point> randPoints; //ds for storing rand points

	std::cout << "Creating random points" << std::endl;
	rand_point_generator(P, chebPoint, num_query_points, var.walk_steps, randPoints, var);
	std::cout << "Created random points" << std::endl;

    std::cout<<"Initialized P..."<<std::endl;

	double naive_time = 0;
	double ann_time = 0;
	int ann_mismatches = 0;
	int ann_explained = 0;
	int not_contained = 0;
	int size = (int)std::ceil(std::pow(1+P.num_of_hyperplanes(), (double)1/2));
	ANNidxArray annIdx;
	ANNdistArray dists;
	annIdx = new ANNidx[size];
	dists = new ANNdist[size];
	ANNpoint queryPt;
	queryPt = annAllocPt(P.dimension());
	double maxMinDist = -1;
	int separatingHyperplaneCount = 0;

	auto it = randPoints.begin();
	for (; it!=randPoints.end(); it++) {
		Point p(P.dimension(), (*it).cartesian_begin(), (*it).cartesian_end());
		auto coords = (*it).cartesian_begin();
		for (int i=0; coords!=(*it).cartesian_end(); i++, coords++) {
			queryPt[i] = (*coords);
		}

		/* Compute min dist to boundary for query point */
		int i=0;
		Point projection = P.project(p, i++);
		double minDist = std::sqrt(((p-CGAL::ORIGIN)-(projection-CGAL::ORIGIN)).squared_length());
		int minDistHyperplane = i;
		for (; i<P.num_of_hyperplanes(); i++) {
			projection = P.project(p, i);
			double tmpDist = std::sqrt(((p-CGAL::ORIGIN)-(projection-CGAL::ORIGIN)).squared_length());
			if (tmpDist<minDist) {
				minDist = tmpDist;
				minDistHyperplane = i;
			}
		}
	
		int naiveSepHyperplane = -1;

		Timer naive_timer;
		bool naive_contains = P.contains_point_naive(p, 0, &naiveSepHyperplane);
		naive_time += naive_timer.elapsed_seconds();

		int nnIndex = -1;
		bool ann_contains;
		Timer ann_timer;
		if (use_jl) {
			ann_contains = P.contains_point_ann_jl(queryPt, annIdx, dists, epsilon, &nnIndex);
		}
		else {
			ann_contains = P.contains_point_ann(queryPt, annIdx, dists, epsilon, &nnIndex);
		}
		ann_time += ann_timer.elapsed_seconds();
	
		if (!naive_contains) {
			not_contained++;
		}
		if (!ann_contains) {
			if (naive_contains) {
				ann_mismatches++;
				if (minDist<epsilon) {
					ann_explained++;
				}
			}
			//Hyperplane naiveHyperplane = P.get_hyperplane(naiveSepHyperplane);
			Hyperplane annHyperplane = P.get_hyperplane(nnIndex);
			//std::cout << "Naive found this separating hyperplane: " << naiveSepHyperplane << " --- " << naiveHyperplane << std::endl;
			std::cout << "Contained naively? " << naive_contains << std::endl;
			std::cout << "ANN found this separating hyperplane: " << nnIndex << " --- " << annHyperplane << std::endl;
			std::cout << "ANN hyperplane has q on negative side? " << annHyperplane.has_on_negative_side(p) << std::endl;
			std::cout << "ANN hyperplane has 0 on positive side? " << annHyperplane.has_on_positive_side(chebPoint) << std::endl;
			std::cout << "Min dist to boundary is " << minDist << " from hyperplane " << minDistHyperplane << " --- " << P.get_hyperplane(minDistHyperplane) << std::endl;
		}
	}
	std::cout << "---Queries inside P--------------------" << std::endl;
	std::cout << "Naive took " << naive_time << "s. " << std::endl;
	std::cout << "ANN took " << ann_time << "s. " << std::endl;
	std::cout << "ANN Mismatch count: " << ann_mismatches << std::endl;
	std::cout << "ANN Explained count: " << ann_explained << std::endl;
	//std::cout << "ANN mismatched point max dist to boundary " << maxMinDist << std::endl;
	//std::cout << "Not contained count: " << not_contained << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	double naive_time2 = 0;
	double ann_time2 = 0;
	double lsh_time2 = 0;
	int lsh_mismatches2 = 0;
	int ann_explained2 = 0;
	int lsh_not_found2 = 0;
	int ann_mismatches2 = 0;
	int not_contained2 = 0;
	it = randPoints.begin();
	for (; it!=randPoints.end(); it++) {
		Point p(P.dimension(), (*it).cartesian_begin(), (*it).cartesian_end());
		auto coords = (*it).cartesian_begin();
		for (int i=0; coords!=(*it).cartesian_end(); i++, coords++) {
			queryPt[i] = (*coords) + epsilon/4;
		}

		/* Compute min dist to boundary for query point */
		int i=0;
		Point projection = P.project(p, i++);
		double minDist = std::sqrt(((p-CGAL::ORIGIN)-(projection-CGAL::ORIGIN)).squared_length());
		int minDistHyperplane = i;
		for (; i<P.num_of_hyperplanes(); i++) {
			projection = P.project(p, i);
			double tmpDist = std::sqrt(((p-CGAL::ORIGIN)-(projection-CGAL::ORIGIN)).squared_length());
			if (tmpDist<minDist) {
				minDist = tmpDist;
				minDistHyperplane = i;
			}
		}
	
		int naiveSepHyperplane = -1;

		Timer naive_timer;
		bool naive_contains = P.contains_point_naive(p, 0, &naiveSepHyperplane);
		naive_time += naive_timer.elapsed_seconds();

		int nnIndex = -1;
		Timer ann_timer;
	
		bool ann_contains;
		if (use_jl) {
			ann_contains = P.contains_point_ann_jl(queryPt, annIdx, dists, epsilon, &nnIndex);
		}
		else {
			ann_contains = P.contains_point_ann(queryPt, annIdx, dists, epsilon, &nnIndex);
		}
		ann_time += ann_timer.elapsed_seconds();
	
		if (naive_contains) {
			not_contained++;
		}
		if (ann_contains) {
			std::cout << "Contained naively? " << naive_contains << std::endl;
			if (!naive_contains) {
				ann_mismatches2++;

				Hyperplane naiveHyperplane = P.get_hyperplane(naiveSepHyperplane);
				std::cout << "Naive found this separating hyperplane: " << naiveSepHyperplane << " --- " << naiveHyperplane << std::endl;
				if (minDist<epsilon) {
					ann_explained2++;
				}
			}
			//Hyperplane annHyperplane = P.get_hyperplane(nnIndex);
			//std::cout << "ANN found this separating hyperplane: " << nnIndex << " --- " << annHyperplane << std::endl;
			//std::cout << "ANN hyperplane has q on negative side? " << annHyperplane.has_on_negative_side(p) << std::endl;
			//std::cout << "ANN hyperplane has 0 on positive side? " << annHyperplane.has_on_positive_side(chebPoint) << std::endl;
			std::cout << "Min dist to boundary is " << minDist << " from hyperplane " << minDistHyperplane << " --- " << P.get_hyperplane(minDistHyperplane) << std::endl;
		}
	}
	delete []annIdx;
	delete []dists;

	std::cout << "---Queries outside P--------------------" << std::endl;
	std::cout << "Naive took " << naive_time2 << "s. " << std::endl;
	std::cout << "ANN took " << ann_time2 << "s. " << std::endl;
	//std::cout << "LSH Mismatch count: " << lsh_mismatches2 << std::endl;
	//std::cout << "ANN not found: " << lsh_not_found2 << std::endl;
	std::cout << "ANN Mismatch count: " << ann_mismatches2 << std::endl;
	std::cout << "ANN Explained count: " << ann_explained2 << std::endl;
	std::cout << "Not contained count: " << not_contained2 << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "---Total--------------------" << std::endl;
	std::cout << "Naive took " << (naive_time+naive_time2) << "s. " << std::endl;
	std::cout << "ANN took " << (ann_time2+ann_time) << "s. " << std::endl;
	//std::cout << "ANN not found: " << (lsh_not_found+lsh_not_found2) << std::endl;
	std::cout << "ANN Mismatch count: " << (ann_mismatches+ann_mismatches2) << std::endl;
	std::cout << "ANN Explained count: " << (ann_explained+ann_explained2) << std::endl;
	std::cout << "Not contained count: " << (not_contained+not_contained2) << std::endl;
	std::cout << "---------------------------------------" << std::endl;
}

int main(const int argc, const char** argv) {
    //Deafault values
    int n, nexp=1, n_threads=1;
    int walk_len;//to be defined after n
    double e=1;
    double exactvol(-1.0);
    bool verbose=false,
         rand_only=false,
         round_only=false,
         file=false,
         round=false,
         NN=false,
         user_walk_len=false,
         linear_extensions=false,
         birk=false,
         rotate=false,
         experiments=true,
         coordinate=true;

    bool membership_test=false;
	bool boundary_test = false;
	double epsilon = 0.1;
	int nqp = 10;
	bool use_jl=true;

    //this is our polytope
    stdHPolytope<double> P;

    if(argc<2) {
        std::cout<<"Use -h for help"<<std::endl;
        exit(-2);
    }

    //parse command line input vars
    for(int i=1; i<argc; ++i) {
        bool correct=false;
        if(!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")) {
            std::cerr<<
                     "Usage:\n"<<
                     "-v, --verbose \n"<<
                     "-rdhr : use random directions HnR, default is coordinate directions HnR\n"
                     "-rand, --rand_only : generates only random points\n"<<
                     "-f1, --file1 [filename type Ax<=b]  [epsilon] [walk length] [threads] [num of experiments]\n"<<
                     //"-f2, --file2 [filename type Ax=b,x>=0] [epsilon] [walk length] [threads] [num of experiments]\n"<<
                     "-fle, --filele : counting linear extensions of a poset\n"<<
                     //"-c, --cube [dimension] [epsilon] [walk length] [threads] [num of experiments]\n"<<
                     "--exact : the exact volume\n"<<
                     "--cube : input polytope is a cube\n"<<
                     "-r, --round : enables rounding of the polytope as a preprocess\n"<<
                     "-ro, --round_only : does only rounding to the polytope\n"<<
                     "-e, --error [epsilon] : the goal error of approximation\n"<<
                     "-w, --walk_len [walk_len] : the random walk length (default 10)\n"<<
                     "-exp [#exps] : number of experiments (default 1)\n"<<
                     "-t, --threads [#threads] : the number of threads to be used\n"<<
                     "-ΝΝ : use Nearest Neighbor search to compute the boundary oracles\n"<<
                     "-birk_sym : use symmetry to compute more random points (only for Birkhoff polytopes)\n"<<
                     std::endl;
            exit(-1);
        }
        if(!strcmp(argv[i],"--cube")) {
            exactvol = std::pow(2,n);
            //exactvol = std::pow(2,n)/std::tgamma(n+1);//factorial of a natural number n is gamma(n+1)
            correct=true;
        }
        if(!strcmp(argv[i],"--exact")) {
            exactvol = atof(argv[++i]);
            correct=true;
        }
        if(!strcmp(argv[i],"-v")||!strcmp(argv[i],"--verbose")) {
            verbose=true;
            std::cout<<"Verbose mode\n";
            correct=true;
        }
        if(!strcmp(argv[i],"-rand")||!strcmp(argv[i],"--rand_only")) {
            rand_only=true;
            std::cout<<"Generate random points only\n";
            correct=true;
        }
        if(!strcmp(argv[i],"-rdhr")) {
            coordinate=false;
            correct=true;
        }
        if(!strcmp(argv[i],"--membership")) {
			std::cout<<"found membership"<<std::endl;
			correct = true;
			boundary_test = false;
			membership_test = true;
        }
        if(!strcmp(argv[i],"--boundary")) {
			std::cout<<"found membership"<<std::endl;
			correct = true;
			boundary_test = true;
			membership_test = false;
        }
        if(!strcmp(argv[i],"--nqp")) {
			correct = true;
            nqp = atoi(argv[++i]);
        }
        if(!strcmp(argv[i],"--epsilon")) {
			correct = true;
            epsilon = atof(argv[++i]);
        }
		if(!strcmp(argv[i],"--ann")) {
			correct = true;
			use_jl = false;
		}
        //reading from file
        if(!strcmp(argv[i],"-f1")||!strcmp(argv[i],"--file1")) {
            file=true;
            std::cout<<"Reading input from file..."<<std::endl;
            std::ifstream inp;
            std::vector<std::vector<double> > Pin;
            inp.open(argv[++i],std::ifstream::in);
            read_pointset(inp,Pin);
            std::cout<<"Read input from file..."<<std::endl;
            //std::cout<<"d="<<Pin[0][1]<<std::endl;
            n = Pin[0][1]-1;
            P.init(Pin);
            std::cout<<"Initialized P..."<<std::endl;
            if (verbose && P.num_of_hyperplanes()<100) {
                std::cout<<"Input polytope: "<<n<<std::endl;
                P.print();
            }
			std::cout <<"Finished reading from file.." << std::endl;
            correct=true;
        }
        /*
            if(!strcmp(argv[i],"-f2")||!strcmp(argv[i],"--file2")){
        			file=true;
              std::ifstream inp;
              std::vector<std::vector<double> > Pin;
              inp.open(argv[++i],std::ifstream::in);
              read_pointset(inp,Pin);
              //std::cout<<"d="<<Pin[0][1]<<std::endl;
              //n = Pin[0][1]-1;
              P.init(Pin);
              P.rref();
              n=P.dimension();
              //if (verbose && P.num_of_hyperplanes()<1000){
        			//	std::cout<<"Input polytope: "<<n<<std::endl;
              //  P.print();
              //}
              correct=true;
            }
        */
        //reading linear extensions and order polytopes
        if(!strcmp(argv[i],"-fle")||!strcmp(argv[i],"--filele")) {
            file=true;
            std::cout<<"Reading input from file..."<<std::endl;
            std::ifstream inp;
            inp.open(argv[++i],std::ifstream::in);
            std::ofstream os ("order_polytope.ine",std::ofstream::out);
            linear_extensions_to_order_polytope(inp,os);

            std::ifstream inp2;
            inp2.open("order_polytope.ine",std::ifstream::in);
            std::vector<std::vector<double> > Pin;
            read_pointset(inp2,Pin);

            //std::cout<<"d="<<Pin[0][1]<<std::endl;
            n = Pin[0][1]-1;
            P.init(Pin);
            //if (verbose && P.num_of_hyperplanes()<100){
            std::cout<<"Input polytope: "<<n<<std::endl;
            //P.print();
            //}
            linear_extensions = true;
            correct=true;
        }
        if(!strcmp(argv[i],"-r")||!strcmp(argv[i],"--round")) {
            round = true;
            correct=true;
        }
        if(!strcmp(argv[i],"-e")||!strcmp(argv[i],"--error")) {
            e = atof(argv[++i]);
            correct=true;
        }
        if(!strcmp(argv[i],"-w")||!strcmp(argv[i],"--walk_len")) {
            walk_len = atof(argv[++i]);
            user_walk_len=true;
            correct=true;
        }
        if(!strcmp(argv[i],"-exp")) {
            nexp = atof(argv[++i]);
            correct=true;
        }
        if(!strcmp(argv[i],"-t")||!strcmp(argv[i],"--threads")) {
            n_threads = atof(argv[++i]);
            correct=true;
        }
        if(!strcmp(argv[i],"-NN")) {
            /*
            if (verbose) std::cout<<"Building search data-srtuctures..."<<std::endl;
            NN=true;
            P.dual(1);
            P.dual(-1);
            */
            std::cout<<"flann software is needed for this option. Experimental feature."
                     <<"Currently under development."<<std::endl;
            correct=true;
        }
        if(!strcmp(argv[i],"-ro")) {
            round_only=true;
            correct=true;
        }
        if(!strcmp(argv[i],"-birk_sym")) {
            birk=true;
            correct=true;
        }
        //rotate the polytope randomly
        if(!strcmp(argv[i],"-rot")) {
            rotate=true;
            correct=true;
        }
        if(correct==false) {
            std::cerr<<"unknown parameters \'"<<argv[i]<<
                     "\', try "<<argv[0]<<" --help"<<std::endl;
            exit(-2);
        }

    }

	if (membership_test) {
		const double err=0.0000000001;
		const double err_opt=0.01;
		//bounds for the cube
		const int lw=0, up=10000, R=up-lw;
		
		/* RANDOM NUMBERS */
		// obtain a time-based seed:
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		// the random engine with this seed
		RNGType rng(seed);
		// standard normal distribution with mean of 0 and standard deviation of 1
		boost::normal_distribution<> rdist(0,1);
		boost::variate_generator< RNGType, boost::normal_distribution<> >
		get_snd_rand(rng, rdist);
		// uniform distribution
		boost::random::uniform_real_distribution<>(urdist);
		boost::random::uniform_real_distribution<> urdist1(-1,1);
    	int rnum = std::pow(e,-2) * 400 * n * std::log(n);
        vars var(rnum,n,walk_len,n_threads,err,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,verbose,rand_only,round,NN,birk,coordinate,use_jl,epsilon);
		membership_main(P,epsilon,nqp,var,use_jl);
	}
	if (boundary_test) {
		const double err=0.0000000001;
		const double err_opt=0.01;
		//bounds for the cube
		const int lw=0, up=10000, R=up-lw;
		
		/* RANDOM NUMBERS */
		// obtain a time-based seed:
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		// the random engine with this seed
		RNGType rng(seed);
		// standard normal distribution with mean of 0 and standard deviation of 1
		boost::normal_distribution<> rdist(0,1);
		boost::variate_generator< RNGType, boost::normal_distribution<> >
		get_snd_rand(rng, rdist);
		// uniform distribution
		boost::random::uniform_real_distribution<>(urdist);
		boost::random::uniform_real_distribution<> urdist1(-1,1);
    	int rnum = std::pow(e,-2) * 400 * n * std::log(n);
        vars var(rnum,n,walk_len,n_threads,err,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,verbose,rand_only,round,NN,birk,coordinate,use_jl,epsilon);
		boundary_main(P,epsilon,nqp,var,use_jl);
	}
	if ((!membership_test) && (!boundary_test)) {
			std::cout<<"Starting old experiments" << std::endl;

    // Set the number of random walk steps
    if(!user_walk_len)
        walk_len=10 + n/10;

    // Timings
    double tstart, tstop;

    /* CONSTANTS */
    //error in hit-and-run bisection of P
    const double err=0.0000000001;
    const double err_opt=0.01;
    //bounds for the cube
    const int lw=0, up=10000, R=up-lw;

    /* RANDOM NUMBERS */
    // obtain a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    // standard normal distribution with mean of 0 and standard deviation of 1
    boost::normal_distribution<> rdist(0,1);
    boost::variate_generator< RNGType, boost::normal_distribution<> >
    get_snd_rand(rng, rdist);
    // uniform distribution
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);

    // If no file specified construct a default polytope
    if(!file) {
        P.init(n);
    }

    // If rotate flag is on rotate the polytope
    if(rotate) {
        rotating(P);
        //P.print();
    }

    // Random walks in K_i := the intersection of the ball i with P
    // the number of random points to be generated in each K_i
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);

    //RUN EXPERIMENTS
    int num_of_exp=nexp;
    double sum=0, sum_time=0;
    double min,max;
    std::vector<double> vs;
    double average, std_dev;
    double Chebtime, sum_Chebtime=double(0);
    NT vol;

    for(int i=0; i<num_of_exp; ++i) {
        std::cout<<"Experiment "<<i+1<<" ";
        //stdHPolytope<double> P_to_test(P);
        tstart = (double)clock()/(double)CLOCKS_PER_SEC;
		std::cout << "Creating ann ds " << std::endl;
		Timer ann_timer;
		Point chebPoint = P.create_point_representation();
		if (use_jl) {
			P.create_ann_jl_ds();
		}
		else {
			P.create_ann_ds();
		}
		double elapsed_ann = ann_timer.elapsed_seconds();
		std::cout << "Created ann ds in " << elapsed_ann << "s " << std::endl;

        // Setup the parameters
        vars var(rnum,n,walk_len,n_threads,err,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,verbose,rand_only,round,NN,birk,coordinate,use_jl,epsilon);


        if(round_only) {
            // Round the polytope and exit
            double round_value = rounding(P,var,var);
            std::cout<<"\n--------------\nRounded polytope\nH-representation\nbegin\n"<<std::endl;
            P.print();
            std::cout<<"end\n--------------\n"<<std::endl;
        } else {
            // Estimate the volume
            vol = volume1_reuse2(P,var,var,Chebtime);
            //if(rotate) vol = std::sqrt(vol);
            //std::cout<<vol<<std::endl;
        }

        double v1 = CGAL::to_double(vol);
        //double v1 = CGAL::to_double(volume1_reuse_test(P_to_test,var,var));
        //double v1 = CGAL::to_double(volume1_reuse_estimete_walk(P_to_test,var,var,Chebtime));

        tstop = (double)clock()/(double)CLOCKS_PER_SEC;
        //double v2 = volume2(P,n,rnum,walk_len,err,rng,get_snd_rand,urdist,urdist1);

        // Statistics
        sum+=v1;
        if(i==0) {
            max=v1;
            min=v1;
        }
        if(v1>max) max=v1;
        if(v1<min) min=v1;
        vs.push_back(v1);
        sum_time +=  tstop-tstart;
        sum_Chebtime += Chebtime;

        //std::cout<<"\t vol= "<<v1<<"\t time= "<<tstop-tstart;
        if(round)
            std::cout<<" (rounding is ON)";
        std::cout<<std::endl;

        //Compute Statistics
        average=sum/(i+1);
        std_dev=0;
        for(std::vector<double>::iterator vit=vs.begin(); vit!=vs.end(); ++vit) {
            std_dev += std::pow(*vit - average,2);
        }
        std_dev = std::sqrt(std_dev/(i+1));

        std::cout.precision(7);

        //MEMORY USAGE
        //struct proc_t usage;
        //look_up_our_self(&usage);

        //Print statistics
        //std::cout<<"\nSTATISTICS:"<<std::endl;
        if (!experiments) {
            std::cout
                    <<"Dimension= "
                    <<n<<" "
                    //<<argv[]<<" "
                    <<"\nNumber of hyperplanes= "
                    <<P.num_of_hyperplanes()<<" "
                    <<"\nNumber of runs= "
                    <<num_of_exp<<" "
                    <<"\nError parameter= "
                    <<e
                    <<"\nTheoretical range of values= "<<" ["
                    <<(1-e)*exactvol<<","
                    <<(1+e)*exactvol<<"] "
                    <<"\nNumber of random points generated in each iteration= "
                    <<rnum<<" "
                    <<"\nRandom walk length= "
                    <<walk_len<<" "
                    <<"\nAverage volume (avg)= "
                    <<average
                    <<"\nmin,max= "
                    " ["
                    <<min<<","
                    <<max<<"] "
                    <<"\nStandard deviation= "
                    <<std_dev<<" "
                    <<"\n(max-min)/avg= "
                    <<(max-min)/average<<" "
                    <<"\nTime(sec)= "
                    <<sum_time/(i+1)<<" "
                    <<"\nTime(sec) Chebyshev= "
                    <<sum_Chebtime/(i+1)<<" "
                    //<<usage.vsize
                    <<std::endl;

            if(exactvol!=-1.0) {
                std::cout
                        <<"\nExact volume= "
                        <<exactvol<<" "
                        <<"\n(vol-avg)/vol= "
                        <<(exactvol-average)/exactvol<<" "
                        <<std::endl;
            }
        } else
            std::cout
                    <<n<<" "
                    //<<argv[]<<" "
                    <<P.num_of_hyperplanes()<<" "
                    <<num_of_exp<<" "
                    <<exactvol<<" "
                    <<e<<" ["
                    <<(1-e)*exactvol<<","
                    <<(1+e)*exactvol<<"] "
                    <<rnum<<" "
                    <<walk_len<<" "
                    <<average<<" ["
                    <<min<<","
                    <<max<<"] "
                    <<std_dev<<" "
                    <<(exactvol-average)/exactvol<<" "
                    <<(max-min)/average<<" "
                    <<sum_time/(i+1)<<" "
                    <<sum_Chebtime/(i+1)<<" "
                    //<<usage.vsize
                    <<std::endl;
    }
    if(linear_extensions)
        std::cout <<"Number of linear extensions= "<<vol*factorial(n)<<std::endl;

    /*
    // EXACT COMPUTATION WITH POLYMAKE
    /*
    std::ofstream polymakefile;
    polymakefile.open("volume.polymake");
    //print_polymake_volfile(C,polymakefile);
    std::cout<<P[0]<<std::endl;
    print_polymake_volfile2(P,polymakefile);
    system ("polymake volume.polymake");
    std::cout<<std::endl;
    */
    //}
}

    return 0;
}


