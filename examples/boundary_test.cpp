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

std::string runTests(stdHPolytope<double>* P, vars& var, int nqp) {
	Point* internalPoint = new Point(P->dimension(), CGAL::ORIGIN);
	json rep;
	Point chebPoint = P->create_point_representation(var, rep, internalPoint);

	CGAL::Random_points_on_sphere_d<Point> rps(P->dimension(), 1);
	// sample nqp points with CDHR

	int totalSteps = 0;
	int maxSteps = -1;
	int minSteps = 101;
	int failed = 0;
	json response;
	response["voronoi"] = rep;
	response["rays"] = json::array();
	for (int i=0; i<nqp; ++i) {
	    std::list<Point> randPoints;
	    rand_point_generator(*P, chebPoint, 1, var.walk_steps, randPoints, var);
		Vector direction = (*rps) - CGAL::ORIGIN;	

		auto it = randPoints.begin();
		Ray r((*it), direction);
		int numberOfSteps = 0;
		bool succeeded = false;
		
		json j;
		P->compute_boundary_intersection(r, &numberOfSteps, &succeeded, 0.001, USE_EXACT, var, j, var.walk_steps, 0);
		if (var.verbose) {
			j["succeeded"] = succeeded;
			response["rays"].push_back(j);
		}
	}
	delete internalPoint;
	return response.dump();
}

int main(int argc, char* argv[]) {
	// parse args
	int n = atoi(argv[1]);
	int d = atoi(argv[2]);
	int diameter = 1000;
	int nqp = 200;
	if (argc>3) {
		diameter = atoi(argv[3]);
	}
	if (argc>4) {
		nqp = atoi(argv[4]);
	}

	// init vars
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	RNGType rng(seed);
	boost::normal_distribution<> rdist(0,1);
	boost::variate_generator< RNGType, boost::normal_distribution<> >
	get_snd_rand(rng, rdist);
	boost::random::uniform_real_distribution<>(urdist);
	boost::random::uniform_real_distribution<> urdist1(-1,1);
   	int rnum = std::pow(0.000001,-2) * 400 * n * std::log(n);
    vars var(rnum,d,d*d*d,1,0.0000001,0,0,0,0,rng,get_snd_rand,
                 urdist,urdist1,true,false,false,false,false,true,0,0.1);

	// create polytope with internal repr
	json results;
	stdHPolytope<double>* P = randomPolytope<double>(n, d, diameter);
	
	results["original"] = runTests(P, var, nqp);

	randomTransformation<stdHPolytope<double> >(P);

	results["transformed"] = runTests(P, var, nqp);

	delete P;
	std::cout << results.dump();

}
