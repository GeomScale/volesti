// VolEsti

// Copyright (c) 2012-2017 Vissarion Fisikopoulos

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

#include "json.hpp"
#include <cstdlib>
#include <vol_rand.h>
#include <rounding.h>
#include <string>
#include <boost/program_options.hpp>


#include <chrono>

namespace funcs2{
class Timer2 {
public:
	Timer2() { start_time = std::chrono::high_resolution_clock::now(); }

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
}

namespace po = boost::program_options; 

json completeTests(stdHPolytope<double>* P, vars& var, int nqp, int k, int l, int algoType, float w) {
	var.verbose = true;
	Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
	json rep;
	json response;
	response["voronoi"] = rep;
	response["rays"] = json::array();
	Point chebPoint = P->create_point_representation(var, rep, internalPoint);
	if (algoType==USE_LSH) {
		P->create_lsh_ds(k, l);
	} 
	else if (algoType==USE_LSHBOX) {
		P->create_lshbox(k, l, w);
	}

	CGAL::Random_points_on_sphere_d<Point> rps(P->dimension(), 1);
	funcs2::Timer2 timer;
	for (int i=0; i<nqp; ++i) {
	    std::list<Point> randPoints;
	    rand_point_generator(*P, chebPoint, 1, var.walk_steps, randPoints, var);
		Vector direction = (*rps) - CGAL::ORIGIN;	

		auto it = randPoints.begin();
		Ray r((*it), direction);
		int numberOfSteps = 0;
		bool succeeded = false;
		json j;

		timer.start();
		auto actualPoint = P->line_intersect(r.source(), r.direction().vector(), false);
		j["actual_time"] = timer.elapsed_seconds();
		j["actual_inside"] = (P->is_in(actualPoint.first)==-1);
		if (!j["actual_inside"]) {
			j["actual_error"] = P->is_in(actualPoint.first);
		}

		timer.start();
		Point appxPoint = P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0.1, algoType, var, j, var.walk_steps, l);
		j["appx_time"] = timer.elapsed_seconds();
		j["appx_path"] = j["steps"];
		j["appx_succeeded"] = succeeded;
		j["appx_steps"] = numberOfSteps;
		j["appx_inside"] = (P->is_in(appxPoint)==-1);
		//if (!j["appx_inside"]) {
			j["appx_error"] = std::sqrt(((actualPoint.first-CGAL::ORIGIN) - (appxPoint-CGAL::ORIGIN)).squared_length());
		//}

		succeeded = false;
		numberOfSteps = 0;
		timer.start();
		Point exactPoint = P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0, USE_EXACT, var, j, var.walk_steps, l);
		j["exact_time"] = timer.elapsed_seconds();
		j["exact_succeeded"] = succeeded;
		j["exact_path"] = j["steps"];
		j["exact_steps"] = numberOfSteps;
		j["exact_inside"] = (P->is_in(exactPoint)==-1);
		if (!j["exact_inside"]) {
			j["exact_error"] = std::sqrt(((actualPoint.first-CGAL::ORIGIN) - (exactPoint-CGAL::ORIGIN)).squared_length());
		}
		j.erase("steps");



		for (auto pit=appxPoint.cartesian_begin(); pit!=appxPoint.cartesian_end(); ++pit) {
			j["appxPoint"].push_back((*pit));
		}
		for (auto pit=exactPoint.cartesian_begin(); pit!=exactPoint.cartesian_end(); ++pit) {
			j["exactPoint"].push_back((*pit));
		}
		for (auto pit=actualPoint.first.cartesian_begin(); pit!=actualPoint.first.cartesian_end(); ++pit) {
			j["actualPoint"].push_back((*pit));
		}
		response["rays"].push_back(j);

		++rps;
	}

	delete internalPoint;
	return response;
}

json simpleTests(stdHPolytope<double>* P, vars& var, int nqp, bool exact=true, int k=-1, int l=-1) {
	Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
	json rep;
	Point chebPoint = P->create_point_representation(var, rep, internalPoint);
	if (!exact) {
		P->create_lsh_ds(k, l);
	}

	CGAL::Random_points_on_sphere_d<Point> rps(P->dimension(), 1);
	// sample nqp points with CDHR

	int totalSteps = 0;
	int maxSteps = -1;
	int minSteps = 101;
	int failed = 0;
	json response;
	response["voronoi"] = rep;
	response["rays"] = json::array();
	auto algo = exact?USE_EXACT:USE_LSH;
	for (int i=0; i<nqp; ++i) {
	    std::list<Point> randPoints;
	    rand_point_generator(*P, chebPoint, 1, var.walk_steps, randPoints, var);
		Vector direction = (*rps) - CGAL::ORIGIN;	

		auto it = randPoints.begin();
		Ray r((*it), direction);
		int numberOfSteps = 0;
		bool succeeded = false;
		
		json j;
		P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0.1, algo, var, j, var.walk_steps, l);
		if (var.verbose) {
			j["succeeded"] = succeeded;
			response["rays"].push_back(j);
		}
	}
	delete internalPoint;
	return response;
}

int main(int argc, char* argv[]) {
	// parse args
	po::options_description desc{"Params"};
	int n;
	int d; 
	int diameter;
	int nqp;
	bool exact = true;
	int k;
	int l;
	int algoType;
	float w;
	desc.add_options()
		("help,h", "Help message")
		(",n", po::value<int>(&n)->required(), "Number of points")
		(",d", po::value<int>(&d)->required(), "Dimension")
		("diameter", po::value<int>(&diameter)->default_value(1000), "Diameter of sphere")
		("nqp", po::value<int>(&nqp)->default_value(100), "Number of query points")
		(",k", po::value<int>(&k)->default_value(-1), "k for LSH")
		(",L", po::value<int>(&l)->default_value(-1), "L for LSH")
		("algo,a", po::value<int>(&algoType)->default_value(USE_LSH), "0 -- FALCONN\n2 -- LSHBOX\n3 -- EXACT")
		(",w", po::value<float>(&w)->default_value(4.0f), "Window of euclidian LSH hash functions")
		("full", "If specified runs every boundary oracle");

	po::variables_map vm; 
	po::store(po::parse_command_line(argc, argv, desc),  vm);
	po::notify(vm);
	//std::cout << "N = " << n << "\nd = " << d << "\ndiam = " << diameter << "\nnqp = " << nqp << "\nk = " << k << "\nL = " << l << std::endl;
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}
	if (k>-1) {
		exact = false;
	}
	// init vars
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	RNGType rng(seed);
	boost::normal_distribution<> rdist(0,1);
	boost::variate_generator< RNGType, boost::normal_distribution<> > get_snd_rand(rng, rdist);
	boost::random::uniform_real_distribution<>(urdist);
	boost::random::uniform_real_distribution<> urdist1(-1,1);
   	int rnum = std::pow(0.000001,-2) * 400 * n * std::log(n);
    vars var(rnum,d,40,1,0.0000001,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,true,false,false,false,false,true,USE_LSHBOX,0.1);

	// create polytope with internal repr
	json results;
	stdHPolytope<double>* P = randomPolytope<double>(n, d, diameter);
	
	if (!vm.count("full")) {
		results["original"] = simpleTests(P, var, nqp, exact, k, l);
	}
	else {
		results["original"] = completeTests(P, var, nqp, k, l, algoType, w);
	}

	randomTransformation<stdHPolytope<double> >(P);	

	if (!vm.count("full")) {
		results["transformed"] = simpleTests(P, var, nqp, exact, k, l);
	}
	else {
		results["transformed"] = completeTests(P, var, nqp, k, l, algoType, w);
	}

	delete P;
	std::cout << results.dump();

}
