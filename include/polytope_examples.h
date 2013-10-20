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


/* Construct a n-CUBE H-REPRESENTATION*/
Polytope cube(int n, NT lw, NT up){	
	Polytope cube;
	std::vector<NT> origin(n,NT(lw));
	for(int i=0; i<n; ++i){
		std::vector<NT> normal;
		for(int j=0; j<n; ++j){
			if(i==j) 
				normal.push_back(NT(1));
			else normal.push_back(NT(0));
		}
		Hyperplane h(Point(n,origin.begin(),origin.end()),
	           Direction(n,normal.begin(),normal.end()));
	  cube.push_back(h);
	}
	std::vector<NT> apex(n,NT(up));
	for(int i=0; i<n; ++i){
		//std::cout<<apex[i]<<" ";
		std::vector<NT> normal;
		for(int j=0; j<n; ++j){
			if(i==j) 
				normal.push_back(NT(-1));
			else normal.push_back(NT(0));
		}
		Hyperplane h(Point(n,apex.begin(),apex.end()),
	           Direction(n,normal.begin(),normal.end()));
	  cube.push_back(h);
	}
	return cube;
}

/* Construct a n-CUBE V-REPRESENTATION*/
V_polytope Vcube(int n, NT lw, NT up){	
	V_polytope cube;
	for(int k=-1; k<2; k+=2){	
		for(int i=0; i<std::pow(2,n-1); ++i){
				//bool bytes[sizeof i];
		    //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
		    //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
		    //          bytes);
				//for(int j=0; j<(sizeof i); ++j)
				boost::dynamic_bitset<> b( n, i );
				std::vector<NT> normal;
				normal.push_back(NT(-1*up*k));
				for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j){
		      if(b[j]) normal.push_back(NT(1*up));
		      else normal.push_back(NT(-1*up));
		    }
		    //Vector normal_v(n,normal.begin(),normal.end());
		    //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
		    //std::cout<<Point(n,normal.begin(),normal.end())<<std::endl;
			  cube.push_back(Point(n,normal.begin(),normal.end()));
		}
	}
	return cube;
}

/* Construct a n-CROSSPOLYTOPE */
Polytope cross(int n, NT lw, NT up){	
	Polytope cross;
	for(int k=-1; k<2; k+=2){
		std::vector<NT> vertex;
		vertex.push_back(NT(k*up));
		for(int j=1; j<n; ++j)
			vertex.push_back(NT(0));
		//std::cout<<Point(n,vertex.begin(),vertex.end())<<std::endl;
		
		for(int i=0; i<std::pow(2,n-1); ++i){
			//bool bytes[sizeof i];
	    //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
	    //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
	    //          bytes);
			//for(int j=0; j<(sizeof i); ++j)
			boost::dynamic_bitset<> b( n, i );
			std::vector<NT> normal;
			normal.push_back(NT(-1*k*up));
			for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j){
	      if(b[j]) normal.push_back(NT(1*up));
	      else normal.push_back(NT(-1*up));
	    }
	    //Vector normal_v(n,normal.begin(),normal.end());
	    //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
			Hyperplane h(Point(n,vertex.begin(),vertex.end()),
			             Direction(n,normal.begin(),normal.end()));
		  cross.push_back(h);
		}
		//std::cout<<"----"<<std::endl;
	}	
	return cross;
}

/* Construct a SKINNY n-CROSSPOLYTOPE */
Polytope cross_skinny(int n, NT lw, NT up){	
	Polytope cross;
	NT sf=pow(2,n);//skinny_factor
	for(int k=-1; k<2; k+=2){
		std::vector<NT> vertex;
		vertex.push_back(NT(k));
		for(int j=1; j<n; ++j)
			vertex.push_back(NT(0));
		//std::cout<<Point(n,vertex.begin(),vertex.end())<<std::endl;
		
		for(int i=0; i<std::pow(2,n-1); ++i){
			//bool bytes[sizeof i];
	    //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
	    //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
	    //          bytes);
			//for(int j=0; j<(sizeof i); ++j)
			boost::dynamic_bitset<> b( n, i );
			std::vector<NT> normal;
			normal.push_back(NT(-1*k*sf));
			for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j){
	      if(b[j]) normal.push_back(NT(1));
	      else normal.push_back(NT(-1));
	    }
	    //Vector normal_v(n,normal.begin(),normal.end());
	    //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
			Hyperplane h(Point(n,vertex.begin(),vertex.end()),
			             Direction(n,normal.begin(),normal.end()));
		  cross.push_back(h);
		}
		//std::cout<<"----"<<std::endl;
	}	
	return cross;
}

//SKINNY 2
Polytope cross_skinny2(int n, NT lw, NT up){	
	Polytope cross;
	NT sf=pow(2,n);//skinny_factor
	for(int k=-1; k<2; k+=2){
		std::vector<NT> vertex;
		vertex.push_back(NT(k*sf));
		for(int j=1; j<n; ++j)
			vertex.push_back(NT(0));
		//std::cout<<Point(n,vertex.begin(),vertex.end())<<std::endl;
		
		for(int i=0; i<std::pow(2,n-1); ++i){
			//bool bytes[sizeof i];
	    //std::copy(static_cast<const bool*>(static_cast<const void*>(&i)),
	    //          static_cast<const bool*>(static_cast<const void*>(&i)) + sizeof i,
	    //          bytes);
			//for(int j=0; j<(sizeof i); ++j)
			boost::dynamic_bitset<> b( n, i );
			std::vector<NT> normal;
			normal.push_back(NT(-1*k/sf));
			for (boost::dynamic_bitset<>::size_type j = 0; j < b.size(); ++j){
	      if(b[j]) normal.push_back(NT(1));
	      else normal.push_back(NT(-1));
	    }
	    //Vector normal_v(n,normal.begin(),normal.end());
	    //std::cout<<Vector(n,normal.begin(),normal.end())<<std::endl;
			Hyperplane h(Point(n,vertex.begin(),vertex.end()),
			             Direction(n,normal.begin(),normal.end()));
		  cross.push_back(h);
		}
		//std::cout<<"----"<<std::endl;
	}	
	return cross;
}

/* Construct a n-CROSS V-REPRESENTATION*/
V_polytope Vcross(int n, NT lw, NT up){	
	V_polytope cross;
	for(int i=0; i<n; ++i){
		std::vector<NT> normal;
		for(int j=0; j<n; ++j){
			if(i==j) 
				normal.push_back(NT(1*up));
			else normal.push_back(NT(0));
		}
	  cross.push_back(Point(n,normal.begin(),normal.end()));
	  //std::cout<<Point(n,normal.begin(),normal.end())<<std::endl;
	}
	for(int i=0; i<n; ++i){
		//std::cout<<apex[i]<<" ";
		std::vector<NT> normal;
		for(int j=0; j<n; ++j){
			if(i==j) 
				normal.push_back(NT(-1*up));
			else normal.push_back(NT(0));
		}
	  //std::cout<<Point(n,normal.begin(),normal.end())<<std::endl;
		cross.push_back(Point(n,normal.begin(),normal.end()));
	}
	return cross;
}


// contruct a n-ball of radius r centered in the origin  
/*
Ball ball(int n, const NT r){
  
  std::vector<Point> P_ball;
  for(int i=0; i<n; ++i){
		std::vector<NT> coords;
		for(int j=0; j<n; ++j){
			if(i==j) 
				coords.push_back(r);
			else coords.push_back(NT(0));
		}
		P_ball.push_back(Point(n,coords.begin(),coords.end()));
	}
	std::vector<NT> extra_coords(n,NT(0));
	extra_coords[0]=NT(-1*r);
	P_ball.push_back(Point(n,extra_coords.begin(),extra_coords.end()));
  Ball B(n,P_ball.begin(),P_ball.end());
	return B;
}
*/

//template <typename T> struct Oracle{
//  sep Sep_Oracle(T &P, Point v);
//};
