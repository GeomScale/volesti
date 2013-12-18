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

int main(const int argc, const char** argv)
{ 
	// VARS
	int n,m;
  if(argc==2){
		n=atoi(argv[1]);
	} else {
		std::cerr<<"Wrong number of arguments"<<std::endl;
    exit(-2);
  }
  m = 2<<(n-1);
  
  std::cout << "cross_"<<n<<".ine\n";
  std::cout << "H-representation\n";
  std::cout << "begin\n";
  std::cout << " " << m << " " << n+1 << " integer\n";
  
  for(int i=0; i<m; ++i){
	//for(int i=m-1; i>=0; --i){
		std::cout << " " << 1 << " ";
		int k=i, j=0;
		while(k!=0){
			int bit = k % 2;
			if(bit == 0)
			  std::cout << "-1 ";
			else
				std::cout << "1 ";
			k = k >> 1;
			++j;
		}
		for(; j<n; ++j)
			std::cout << -1 << " ";
		std::cout << "\n";
	}
	
	std::cout << "end\ninput_incidence" << std::endl;	
    
	return 0;  	
}
