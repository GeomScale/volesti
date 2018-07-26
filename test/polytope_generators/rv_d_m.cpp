// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2017 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

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
  
  std::cout << "cube_"<<n<<".ine\n";
  std::cout << "H-representation\n";
  std::cout << "begin\n";
  std::cout << " " << 2*n << " " << n+1 << " integer\n";
  
  for(int i=0; i<n; ++i){
		std::cout << " " << 1 << " ";
		for(int j=0; j<n; ++j){
			if(i==j) 
			  std::cout << 1 << " ";
			else std::cout << 0 << " ";
		}
		std::cout << "\n";
	}
	for(int i=0; i<n; ++i){
		std::cout << " " << 1 << " ";
		for(int j=0; j<n; ++j){
			if(i==j) 
			  std::cout << -1 << " ";
			else std::cout << 0 << " ";
		}
		std::cout << "\n";
	}
	std::cout << "end\ninput_incidence" << std::endl;	
    
	return 0;  	
}
