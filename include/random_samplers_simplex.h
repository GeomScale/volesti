// VolEsti

// Copyright (c) 2017 Tolis Chalkis

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

/*

Sampling Uniformly from Unit or an Arbitrary Simplex given in V-representation.

By Tolis Chalkis, N.K. Univ. Athens, Greece, 2017.

-----For sampling from Unit Simplex we give two implemetations-----

Both of them are based on Noah A.Smith and Roy W. Tromble "Sampling Uniformly from the Unit Simplex".

For SampleUnitSimplex01() we implement a variation of Filter Bloom for dimensions greater than 60 to
ensure a distinct choice of integers. For dimensions lower than 60 we use std::vector and a linear
search for each integer's uniqness in every step, because is faster than the Filter.

For SampleUnitSimplex02() we implement a hash function (using division method) for all dimensions to
ensure a distinct choice of integers. Using std::vector to ensure a distinct choice is faster for d<130.

std::unordered_map or std::set are much more slower to use in order to ensure a distinct choice of integers.
But we give the implemetations in SampleUnitSimplex03() and SampleUnitSimplex04() as comments at the end.
Using std::vector to ensure a distinct choice is faster than std::unordered_map for d<400 and than
std::set for d<600.

SampleUnitSimplex01() is faster than all the others for every dimension.


-----For sampling from an arbitrary Simplex we give three implemetations--------

SampleArbSimplex01() is based on Art Owen "Monte Carlo theory, methods and examples", section 5.4 and
is the fastest implematation. To ensure a distinct choice of integers we use the same variation of 
Filter Bloom, as before, for dimensions greater than 60. For dimensions lower than 60 we use std::vector
and a linear search for each integer's uniqness in every step.

SampleArbSimplex02() is sampling from Unit Simplex and then applying a uniformly-preserving linear
transformation on these points mapping them to the arbitrary Simplex. SampleUnitSimplex02() is faster
than SampleUnitSimplex03().

SampleArbSimplex03() is based on Christian Grimme "Picking a Uniformly Random Point from an Arbitrary
Simplex". It's quite different algorithm than the other two but is the slowest implemetation.


-----All the implementations check only if the number of vertices is equal to d+1. Then the algorithms
assume that the polytope is a Simplex.


*****************************************************************
-------------------How to call the functions---------------------
*****************************************************************


------------ Sampling from Unit Simplex------------

The functions take three Arguments.

The first Argument is an integer indicates the dimension.
The second Argument is an integer indicates the number of random points requested.
The third Argument should be an empty Point_d-vector. All the random points will be stored in it.

For Example, if:

- std::vector<Point_d> Points; Is the empty Point_d-vector

Then the call --> SampleUnitSimplex01(3,5000,Points); <-- Would generate 5000 uniformly random points
in the Unit Simplex for d=3 and store them to the vector "Points".



-------- Sampling from an Arbitrary Simplex-------

The functions take four arguments. First we need to have the vertices in a vector that contains Points
in d dimension.

Then the first two Arguments are the iterators pointing at the beginning and the end of the vertices-vector.
The third Argument is an integer indicates the number of random points requested.
The fourth Argument should be an empty Point_d-vector. All the random points will be stored in it.

For Example, if:

- std::vector<Point_d> v; Is the vector that contains the vertices
- std::vector<Point_d> Points; Is the empty Point_d-vector

Then the call --> SampleArbSimplex01(v.begin(),v.end(),5000,Points); <-- Would generate 5000 uniformly
random points and store them to the vector "Points".


--------- For more details and examples on how to call the functions see the sample_simplex_main.cpp

*/



#include <CGAL/Cartesian_d.h>
#include <vector>
#include "boost/random.hpp"

//#include <tr1/unordered_map>
//#include <set>


typedef boost::mt19937 			RNGType;
typedef CGAL::Cartesian_d<double> 	        Kernel; 
typedef Kernel::Point_d 		Point_d;




void SampleUnitSimplex01(int dim, int num, std::vector<Point_d> &points){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<bool> filter;
	
	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	x_vec.assign(dim+1,0);

	if (dim>60){
		bool t1=true,t2=true;
		pr=3*dim+1;
		if(pr%2==0) pr+=1;

		while(t1){
			t2=true;
			divisors=(int)floor(sqrt((double)pr))+1;
			for (i=3; i<divisors+1; i+=2){
				if (pr%i==0){
					t2=false;
					break;
				}
			}
			if(t2) break;
			pr+=2;
		}
		
		// Generate the number of points requested
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and the filter equal to true
			y.assign(dim,0); filter.assign(pr,true);
			pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){
				x_rand = uidist(rng);// md=x_rand%pr;

				// Check if this integer is the first that is mapped to the specific filter's position
				if ( filter[x_rand%pr] ){
					filter[x_rand%pr]=false; pointer++;
					x_vec[pointer]=x_rand;
				}
			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );
	
			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}else{
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and clear the integers' list
			y.assign(dim,0); pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){

				x_rand = uidist(rng);
				// Check if this integer is selected first time
				if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand) == x_vec.begin()+pointer ){
					pointer++;
					x_vec[pointer]=x_rand;
				}

			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );

			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}
	return;
}








void SampleUnitSimplex02(int dim, int num, std::vector<Point_d> &points){
	
	int j,i,k,x_rand,M=2147483647,pr,divisors,index;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<bool> filter;
	
	bool t1=true,t2=true;

	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	// Find the nearest prime>d/3
	pr=(int)floor(dim/3)+1;
	if(pr%2==0) pr+=1;

	while(t1){
		t2=true;
		divisors=(int)std::floor(std::sqrt((double)pr))+1;
		for (i=3; i<divisors+1; i+=2){
			if (pr%i==0){
				t2=false;
				break;
			}
		}
		if(t2) break;
		pr+=2;
	}
		
	std::vector< std::vector<int> > table(pr);
	x_vec.assign(dim+1,0);

	// Generate the number of points requested
	for (i=0; i<num; i++){

		// Set all the point's coordinates equal to zero and the filter equal to true
		y.assign(dim,0); index=0;
	
		for (j=0; j<pr; j++){
			table[j].clear();
		}
		
		// Generate d distinct integers
		while ( index<dim ){
			x_rand = uidist(rng);
			// Check if this integer is selected first time
			if ( std::find(table[x_rand%pr].begin(), table[x_rand%pr].end(), x_rand) == table[x_rand%pr].end() ){
				table[x_rand%pr].push_back(x_rand); index++;
				x_vec[index]=x_rand;
			}
		}

		// Sort the integers' list
		std::sort( x_vec.begin(), x_vec.end() );
		
		// Construct the point's coordinates
		for(j=0; j<dim; j++){
			y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
		}

		// Define the new point and add it to the point-vector
		points.push_back(Point_d(dim,y.begin(),y.end()));
		
	}
	return;
}








void SampleArbSimplex01(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end, int num, std::vector<Point_d> &points){
	
	int n=std::distance(it_beg,it_end),dim,j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<bool> filter;
	double Xj;

	Point_d p0=*it_beg;
	if (n>=2){   // Check if points could define a simplex
		if (p0.dimension()!=n-1){   // Check if the polytope is a Simplex
			std::cout<<"This Polytope is not a Simplex. The problem is open."<<std::endl;
			return;
		}else{
			dim=p0.dimension();
		}
	}else{
		std::cout<<"We need more points"<<std::endl;
		return;
	}

	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	x_vec.assign(dim+2,0);  x_vec[dim+1]=M;
	
	if (dim>60){
		bool t1=true,t2=true;

		// Find the nearest prime>3d
		pr=3*dim+1;
		if(pr%2==0) pr+=1;

		while(t1){
			t2=true;
			divisors=(int)std::floor(std::sqrt((double)pr))+1;
			for (i=3; i<divisors+1; i+=2){
				if (pr%i==0){
					t2=false;
					break;
				}
			}
			if(t2) break;
			pr+=2;
		}
		
		// Generate the number of points requested
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and the filter equal to true
			y.assign(dim,0); filter.assign(pr,true);
			pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){
				x_rand = uidist(rng);

				// Check if this integer is the first that is mapped to the specific filter's position
				if ( filter[x_rand%pr] ){
					filter[x_rand%pr]=false; pointer++;
					x_vec[pointer]=x_rand;
				}
			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );
		
			// Construct the point's coordinates
			for(j=0; j<dim+1; j++){
				Point_d pk=*(it_beg+j);
				Xj=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for (k=0; k<dim; k++){
					y[k]+=Xj*pk[k];
				}
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		}
		
	}else{
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and clear the integers' list
			y.assign(dim,0); pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){

				x_rand = uidist(rng);
				// Check if this integer is selected first time
				if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand) == x_vec.begin()+pointer ){
					pointer++;
					x_vec[pointer]=x_rand;
				}

			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );

			// Construct the point's coordinates
			for(j=0; j<dim+1; j++){
				Point_d pk=*(it_beg+j);
				Xj=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for (k=0; k<dim; k++){
					y[k]+=Xj*pk[k];
				}
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}
	return;
}









void SampleArbSimplex02(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end, int num, std::vector<Point_d> &points){
	
	int n=std::distance(it_beg,it_end),dim,j,i,k,x_rand,M=2147483647,pr,divisors,pointer;  // M is the largest possible integer
	std::vector<int> x_vec;
	std::vector<double> y;
	std::vector<bool> filter;
	double coord_i;
	
	Point_d p0=*it_beg;
	if (n>=2){   // Check if points could define a simplex
		if (p0.dimension()!=n-1){   // Check if the polytope is a Simplex
			std::cout<<"This Polytope is not a Simplex. The problem is open."<<std::endl;
			return;
		}else{
			dim=p0.dimension();
		}
	}else{
		std::cout<<"We need more points"<<std::endl;
		return;
	}

	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	x_vec.assign(dim+1,0);

	
	if (dim>60){

		bool t1=true,t2=true;

		// Find the nearest prime>3d
		pr=3*dim+1;
		if(pr%2==0) pr+=1;

		while(t1){
			t2=true;
			divisors=(int)std::floor(std::sqrt((double)pr))+1;
			for (i=3; i<divisors+1; i+=2){
				if (pr%i==0){
					t2=false;
					break;
				}
			}
			if(t2) break;
			pr+=2;
		}
		
		// Generate the number of points requested
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and the filter equal to true
			y.assign(dim,0); filter.assign(pr,true);
			pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){
				x_rand = uidist(rng);// md=x_rand%pr;

				// Check if this integer is the first that is mapped to the specific filter's position
				if ( filter[x_rand%pr] ){
					filter[x_rand%pr]=false; pointer++;
					x_vec[pointer]=x_rand;
				}
			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );
		
			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				Point_d pk=*(it_beg+j+1);
				coord_i=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for(k=0; k<dim; k++){
					y[k]+=coord_i*(pk[k]-p0[k]);
				}
				y[j]+=p0[j];
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		}
		
	}else{
		for (i=0; i<num; i++){

			// Set all the point's coordinates equal to zero and clear the integers' list
			y.assign(dim,0); pointer=0;
		
			// Generate d distinct integers
			while ( pointer<dim ){

				x_rand = uidist(rng);
				// Check if this integer is selected first time
				if ( std::find(x_vec.begin(), x_vec.begin()+pointer, x_rand) == x_vec.begin()+pointer ){
					pointer++;
					x_vec[pointer]=x_rand;
				}

			}

			// Sort the integers' list
			std::sort( x_vec.begin(), x_vec.end() );

			// Construct the point's coordinates
			for(j=0; j<dim; j++){
				Point_d pk=*(it_beg+j+1);
				coord_i=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
				for(k=0; k<dim; k++){
					y[k]+=coord_i*(pk[k]-p0[k]);
				}
				y[j]+=p0[j];
			}

			// Define the new point and add it to the point-vector
			points.push_back(Point_d(dim,y.begin(),y.end()));
		
		}
	}
	return;
}










void SampleArbSimplex03(std::vector<Point_d>::iterator it_beg, std::vector<Point_d>::iterator it_end, int num, std::vector<Point_d> &points){
	
	int n=std::distance(it_beg,it_end),j,dim,i,k,m,x_rand,M=2147483647;  // M is the largest possible integer
	std::vector<double> y;
	std::vector<double> l;
	double bi;
	Point_d p0=*it_beg;
	if (n>=2){   // Check if points could define a simplex
		if (p0.dimension()!=n-1){   // Check if the polytope is a Simplex
			std::cout<<"This Polytope is not a Simplex. The problem is open."<<std::endl;
			return;
		}else{
			dim=p0.dimension();
		}
	}else{
		std::cout<<"We need more points"<<std::endl;
		return;
	}

	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(0,M);
	l.assign(dim+2,0); l[0]=1;

	// Generate the number of points requested
	for (m=0; m<num; m++){    

		y.assign(dim,0);    // Set all the point's coordinates equal to zero

		//Construct the weights
		for (j=1; j<n; j++){
			x_rand = uidist(rng);
			l[j]=std::pow (((double)x_rand)/((double)M), 1/((double)(n-j)) );  
		}
	
		for (i=1; i<dim+2; i++){
			// Construct the coefficients of the vertices
			bi=1-l[i];
			for (j=0; j<i; j++){
				bi*=l[j];
			}

			// Construct the coordinates of the point
			Point_d pi=*(it_beg+i-1);
			for (k=0; k<dim; k++){
				y[k]+=bi*pi[k];
			}
		}

		// Define the new point and add it to the point-vector
		points.push_back(Point_d(dim,y.begin(),y.end()));

	}
	
	return;
}





//---------------------------------------------------------------------------------------//
/*    The two slower implementations for the Unit Simplex are given below



void SampleUnitSimplex03(int dim, int num, std::vector<Point_d> &points){
	
	int j,i,k,x_rand,M=2147483647,pointer;
	std::tr1::unordered_map<int, bool>	h_table;
	std::tr1::unordered_map<int, bool>::iterator got;
	std::vector<int> x_vec;
	std::vector<double> y;
	

	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	x_vec.assign(dim+1,0);

	// Generate the number of points requested
	for (i=0; i<num; i++){

		// Set all the point's coordinates equal to zero
		y.assign(dim,0); h_table.clear();
		pointer=0;
		
		// Generate d distinct integers
		while ( pointer<dim ){
			x_rand = uidist(rng);
			got = h_table.find(x_rand);

			if ( got == h_table.end() ){   // Checks if the integer is chosen for first time
				h_table.insert(std::pair<int,bool> (x_rand,true)); pointer++;
				x_vec[pointer]=x_rand;
			}
		}

		// Sort the integers' list
		std::sort( x_vec.begin(), x_vec.end() );
		
		// Construct the point's coordinates
		for(j=0; j<dim; j++){
			y[j]=((double)(x_vec[j+1]-x_vec[j])) / ((double)M);
		}

		// Define the new point and add it to the point-vector
		points.push_back(Point_d(dim,y.begin(),y.end()));
		
	}
	return;
}




void SampleUnitSimplex04(int dim, int num, std::vector<Point_d> &points){
	
	int j,i,k,x_rand,M=2147483647;   // M is the largest possible integer
	std::set<int> x_vec;
	std::set<int>::iterator it;
	std::set<int>::iterator it2;
	std::vector<double> y;

	RNGType rng;
	boost::random::uniform_int_distribution<> uidist(1,M);

	// Generate the number of points requested
	for (i=0; i<num; i++){

		// Set all the point's coordinates equal to zero and clear the integers' set
		y.assign(dim,0);
		x_vec.clear(); x_vec.insert(0);

		// Generate d distinct integers
		while (x_vec.size()<(dim+1)){
			x_rand = uidist(rng);
			x_vec.insert(x_rand);  // The integer is inserted only if is selected fist time
		}

		it2=x_vec.begin(); it2++;
		it=x_vec.begin();

		// Construct the point's coordinates
		for(j=0; j<dim; j++){
			y[j]=((double)(*it2-*it)) / ((double)M);
			it2++; it++;
		}

		// Define the new point and add it to the point-vector
		points.push_back(Point_d(dim,y.begin(),y.end()));
		
	}
	return;
}





*/
