// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2017 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

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
