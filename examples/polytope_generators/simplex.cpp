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

#include <iostream>

int main(const int argc, const char** argv)
{ 
	// VARS
	int n;
  if(argc==2){
		n=atoi(argv[1]);
	} else {
		std::cerr<<"Wrong number of arguments"<<std::endl;
    exit(-2);
  }
  
  std::cout << "simplex_"<<n<<".ine\n";
  std::cout << "H-representation\n";
  std::cout << "begin\n";
  std::cout << " " << n+1 << " " << n+1 << " integer\n";
  
  for(int i=0; i<n; ++i){
		std::cout << 0 << " ";
		for(int j=0; j<n; ++j){
			if(i==j) 
			  std::cout << 1 << " ";
			else std::cout << 0 << " ";
		}
		std::cout << "\n";
	}
	std::cout << 1 << " ";
	for(int j=0; j<n; ++j){
		  std::cout << -1 << " ";
	}
	std::cout << "\n";
	std::cout << "end\ninput_incidence" << std::endl;	
    
	return 0;  	
}
