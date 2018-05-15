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

###################################################
--------------------- Main ------------------------
###################################################

This main file could help in execution for fuctions in random_samplers_simplex.h 

------------ Sampling from Unit Simplex------------

The command:
./executable unit1 12 30000
Will call SampleUnitSimplex01() for d=12 and 30000 random points requested.

The command:
./executable unit2 15 20000
Will call SampleUnitSimplex02() for d=12 and 20000 random points requested.


-------- Sampling from an Arbitrary Simplex-------

For sampling from an arbitrary Simplex we have to give the vertices through a txt file.
At the first line we should only declare the dimension d.
In d+1 lines below we should declare the Points-vertices of the Simplex. For example the txt:

***************
* 3           *
* 13 4 -27    *
* 12 57 127   *
* 12.3 -78 91 *
* -7 89 21    *
***************

Is a valid txt for this main. 

The command:
./executable arb1 vertices.txt 100000
Will call SampleArbSimplex01() for 100000 random points requested.

The command:
./executable arb2 vertices.txt 70000
Will call SampleArbSimplex02() for 70000 random points requested.

The command:
./executable arb3 vertices.txt 6000
Will call SampleArbSimplex03() for 6000 random points requested.


------- We create a txt file that contains the random points (random_points.txt) --------

The command -------> ./executable -help <------- Will print all the options.

*/


#include <fstream>
#include <random_samplers_simplex.h>


typedef double NT;
typedef CGAL::Cartesian_d<NT> 	      Kernel; 
typedef Kernel::Point_d Point_d;


int main(int argc, char* argv[]){

	std::vector<NT> newp;
	std::vector<Point_d> v;
	std::vector<Point_d> Points;
	std::vector<Point_d>::iterator iter;
	int dim,num,i;

	if ( argc!=4 && argc!=2 ){
		printf ("%d\n",argc);
		std::cout<<"Wrong inputs! Use -help for help"<<std::endl;

	}else if( !strcmp(argv[1],"-help") ){

          std::cerr<<
		"Usage:\n"<<
		"-arb1 : Sampling uniformly from an arbitrary Simplex, using an Owen's algorithm Implementation\nCommand : ./exe -arb1 vertices.txt number_of_points\n\n"<<
		"-arb2 : Sampling uniformly from an arbitrary Simplex, applying a linear transformation on Smith and Tromble's algorithm's random points in Unit Simplex\nCommand : ./exe -arb2 vertices.txt number_of_points\n\n"<<
		"-arb3 : Sampling uniformly from an arbitrary Simplex, using a Grimme's algorithm Implementation\nCommand : ./exe -arb3 vertices.txt number_of_points\n\n"<<
		"-unit1 : Sampling Uniformly from the Unit Simplex with Smith and Tromble's algorithm, using a Bloom filter variation to ensure a distinct choice of integers\nCommand : ./exe -unit1 dimension number_of_points\n\n"<<
		"-unit2 : Sampling Uniformly from the Unit Simplex with Smith and Tromble's algorithm, using a hash function to ensure a distinct choice of integers \nCommand : ./exe -unit2 dimension number_of_points\n"<<
	  std::endl;
          exit(-1);

	}else if( !strcmp(argv[1],"-arb1") || !strcmp(argv[1],"-arb2") || !strcmp(argv[1],"-arb3") ){

		num=atoi(argv[3]);

		std::ifstream inputFile;
		std::string inputFileName = argv[2]; //reading the input txt file
		inputFile.open(inputFileName.c_str()); //converting the input to a c-style string
	
		std::istream_iterator< NT >  input_begin( inputFile );
		std::istream_iterator< NT >  input_end;
		for(std::istream_iterator< NT > it = input_begin; it != input_end; ++it){
			newp.push_back(*it);
		}
		dim=newp[0];    // Get the dimension
		int n=(newp.size()-1)/dim;  // Get the number of the points
		v.reserve (n);
		for(i=0;i<n;i++){
			v.push_back(Point_d(dim,newp.begin()+i*dim+1,newp.begin()+i*dim+dim+1));  // Insert points to the vector
		}
		inputFile.close();
	
		if( !strcmp(argv[1],"-arb1") ){
			SampleArbSimplex01(v.begin(),v.end(),num,Points);
		}else if( !strcmp(argv[1],"-arb2") ){
			SampleArbSimplex02(v.begin(),v.end(),num,Points);
		}else{
			SampleArbSimplex03(v.begin(),v.end(),num,Points);
		}

		std::ofstream outputFile;
		outputFile.open("random_points.txt");   // Create the txt file that contains the random points

		for (iter=Points.begin(); iter!=Points.end(); iter++){
			outputFile<< *iter <<std::endl;
		}

		outputFile.close();

	}else if( !strcmp(argv[1],"-unit1") || !strcmp(argv[1],"-unit2") ){

		num=atoi(argv[3]); dim=atoi(argv[2]);

		if ( !strcmp(argv[1],"-unit1") ){
			SampleUnitSimplex01(dim,num,Points);
		}else{
			SampleUnitSimplex02(dim,num,Points);
		}
		
		std::ofstream outputFile;
		outputFile.open("random_points.txt");    // Create the txt file that contains the random points

		for (iter=Points.begin(); iter!=Points.end(); iter++){
			outputFile<< *iter <<std::endl;
		}

		outputFile.close();

	}else{

		std::cout<<"Wrong inputs! Use -help for help"<<std::endl;

	}

	return 0;
	
}
