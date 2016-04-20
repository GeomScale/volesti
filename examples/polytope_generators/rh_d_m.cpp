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

#include <iostream>
#include <time.h>

int main(const int argc, const char** argv)
{ 
	// VARS
	int d, m;
  if(argc==3){
		d=atoi(argv[1]);
		m=atoi(argv[2]);
	} else {
		std::cerr<<"Wrong number of arguments"<<std::endl;
    exit(-2);
  }
  
  std::cout << "rh_"<<d<<"_"<<m<<".ine\n";
  std::cout << "H-representation\n";
  std::cout << "begin\n";
  std::cout << " " << m << " " << d+1 << " integer\n";
  
  srand (time(NULL));
  
  for(int i=1; i<m; ++i){
		std::cout << " " << 1000 << " ";
		for(int j=0; j<d; ++j){
			if(std::rand()%2==1)
				std::cout << "-";
			std::cout << std::rand()%1000 << " ";
		}
		std::cout << "\n";
	}
	std::cout << "end\ninput_incidence" << std::endl;	
    
	return 0;  	
}
