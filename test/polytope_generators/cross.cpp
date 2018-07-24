// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2017 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

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
