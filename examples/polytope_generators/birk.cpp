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
// Public License.  If you did not receive this file along with HeaDDaCHe,
// see <http://www.gnu.org/licenses/>.
// 
// Developer: Vissarion Fisikopoulos

#include <iostream>
#include <math.h>

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
  int m = pow(n,2);
  int d = pow(n-1,2)+1;
  
  std::cout << "birk_"<<n<<".ine\n";
  std::cout << "H-representation\n";
  std::cout << "begin\n";
  std::cout << " " << m << " " << d << " integer\n";
    
  std::cout << -1*(n-2) << " ";  
  for(int j=1; j< d; ++j)
    std::cout << "1 ";
  std::cout << "\n";
    
  for(int i=0; i<n-1; ++i){
		std::cout << "1 ";
		for(int j=1; j< d; ++j){
			if(j%(n-1) == i){
		    std::cout << "-1 ";
		  }
		  else std::cout << " 0 ";
		}
		std::cout << "\n";
	}
	
	for(int i=0; i<n-1; ++i){
		std::cout << "1 ";
		for(int j=0; j< d-1; ++j){
			if(j/(n-1) == i){
		    std::cout << "-1 ";
		  }
		  else std::cout << " 0 ";
		}
		std::cout << "\n";
	}
	
	for(int i=0; i<d-1; ++i){
		std::cout << "0 ";
		for(int j=0; j< d-1; ++j){
			if(j == i){
		    std::cout << " 1 ";
		  }
		  else std::cout << " 0 ";
		}
		std::cout << "\n";
	}
	std::cout << "end\ninput_incidence" << std::endl;	
   
 
   
	return 0;  	
}
