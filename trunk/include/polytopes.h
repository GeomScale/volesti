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

#ifndef POLYTOPES_H
#define POLYTOPES_H


// my V polytope class
template <typename K>
class stdHPolytope{
	private:
	  typedef std::vector<K>        stdCoeffs;
	  typedef std::vector<stdCoeffs>  stdMatrix;
	
	public: 
	  // default constructor: cube(d)
	  stdHPolytope(int d): _d(d) {
			for(int i=0; i<d; ++i){
				stdCoeffs coeffs;
				coeffs.push_back(K(1));
				for(int j=0; j<d; ++j){
					if(i==j) 
					  coeffs.push_back(K(1));
					else coeffs.push_back(K(0));
				}
			_A.push_back(coeffs);
		  }
		  for(int i=0; i<d; ++i){
				stdCoeffs coeffs;
				coeffs.push_back(K(1));				
				for(int j=0; j<d; ++j){
					if(i==j) 
					  coeffs.push_back(K(-1));
					else coeffs.push_back(K(0));
				}
			_A.push_back(coeffs);
		  }
	  }
	  
	  int print() {
			for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
				for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit)
				  std::cout<<*lit<<" ";
		    std::cout<<std::endl;
		  }			
			return 0;
		}
	  
	  int size() {
			return _A.size();
		}
	  
	  /*
	  int is_in(stdCoeffs p) {
			for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
				typename stdCoeffs::iterator lit,pit;
				pit=p.begin(); 
				lit=mit->coefficients_begin();
				K sum=(*lit);
				++lit;
				for( ; lit<mit->coefficients_end() ; ++lit){
					//std::cout << *lit << " " << *pit <<std::endl;
					sum += *lit * (*pit);
				}
				//std::cout<<sum<<std::endl;
				if(sum<K(0))
					return mit-_A.begin();
			}
			return -1;
		}
	  */
	  
	  int is_in(Point p) {
			for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
				typename stdCoeffs::iterator lit;
				Point::Cartesian_const_iterator pit;
				pit=p.cartesian_begin(); 
				lit=mit->begin();
				K sum=(*lit);
				++lit;
				for( ; lit<mit->end() ; ++lit, ++pit){
					//std::cout << *lit << " " << *pit <<std::endl;
					sum += *lit * (*pit);
				}
				
				//std::cout<<sum<<std::endl;
				if(sum<K(0))
					return mit-_A.begin();
			}
			return -1;
		}
	  
	  // compute intersection point of ray starting from r and pointing to v
	  // with polytope discribed by _A
	  std::pair<Point,Point> line_intersect(Point r, 
                         Vector v){
		  //std::cout<<"line-polytope-intersection"<<std::endl;
		  K lamda=0;
		  K min_plus=0, max_minus=0;
		  bool min_plus_not_set=true;
		  bool max_minus_not_set=true;
		  for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
				typename stdCoeffs::iterator cit;
				Point::Cartesian_const_iterator rit;
				rit=r.cartesian_begin(); 
				Point::Cartesian_const_iterator vit;
				vit=v.cartesian_begin();
				cit=ait->begin();
				K sum_nom=(*cit);
				++cit;
				K sum_denom=K(0); 
				//std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
				//         v.cartesian_begin()-v.cartesian_end()<<std::endl;
				for( ; cit < ait->end() ; ++cit, ++rit, ++vit){
					//std::cout << sum_nom << " " << sum_denom <<std::endl;
					//std::cout << int(rit-r.cartesian_begin()) << " " << int(vit-v.cartesian_begin()) <<std::endl;
					sum_nom -= *cit * (*rit);
					sum_denom += *cit * (*vit);
				}
			  if(sum_denom==K(0)){
			    std::cout<<"div0"<<std::endl;
			  }
			  else{
			    lamda = sum_nom/sum_denom;
				  if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
				  if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
				  if(lamda<min_plus && lamda>0) min_plus=lamda;
				  if(lamda>max_minus && lamda<0) max_minus=lamda;
				}			  
				//std::cout<<r+(lamda*v)<<"\n"<<lamda<<std::endl;
			}
			/*
			std::cout<<"lmin,lmax= "<<max_minus<<" "<<min_plus<<std::endl;
			std::cout<<"r= "<<r<<std::endl;
			std::cout<<"v= "<<v<<std::endl;
			std::cout<<"p1= "<<r+(min_plus*v)<<std::endl;
			std::cout<<"p2= "<<r+(max_minus*v)<<std::endl;
			*/
			return std::pair<Point,Point> (r+(min_plus*v),r+(max_minus*v));
		}
	  
	  std::pair<NT,NT> line_intersect_coord(Point &r,
														              int rand_coord){
			//std::cout<<"line-polytope-intersection"<<std::endl;
		  K lamda=0;
		  //std::vector<NT> new_lamdas(_A.size());
		  //std::vector<NT> new_lamdas;
		  K min_plus=0, max_minus=0;
		  bool min_plus_not_set=true;
		  bool max_minus_not_set=true;
		  //std::vector<NT>::iterator lamdait = lamdas.begin();
		  for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
				typename stdCoeffs::iterator cit;
				Point::Cartesian_const_iterator rit;
				rit=r.cartesian_begin(); 
				//Point::Cartesian_const_iterator vit;
				//vit=v.cartesian_begin();
				cit=ait->begin();
				K sum_nom=(*cit);
				++cit;
				//here we just need the "rand_coord" coordinate of c
				//std::cout<<*(cit+rand_coord)<<"c= "<<std::endl;
				//for(typename stdCoeffs::iterator cit2=ait->begin() ; cit2 < ait->end() ; ++cit2){
				//  std::cout<<*cit2<<" ";
				//}
				//std::cout<<std::endl;
				K sum_denom= *(cit+rand_coord); 
				//std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
				//         v.cartesian_begin()-v.cartesian_end()<<std::endl;
				for( ; cit < ait->end() ; ++cit, ++rit){
					//std::cout << sum_nom << " " << sum_denom <<std::endl;
					//std::cout << int(rit-r.cartesian_begin()) << " " << int(vit-v.cartesian_begin()) <<std::endl;
					sum_nom -= *cit * (*rit);
					//sum_denom += *cit * (*vit);
				}
				//std::cout << sum_nom << " / "<< sum_denom<<std::endl;
			  if(sum_denom==K(0)){
			    //std::cout<<"div0"<<std::endl;
			    ;
			  }
			  else{
			    lamda = sum_nom*(1/sum_denom);
			    //lamdas[ait-_A.begin()] = lamda;
			    
				  if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
				  if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
				  if(lamda<min_plus && lamda>0) min_plus=lamda;
				  if(lamda>max_minus && lamda<0) max_minus=lamda;
				}			  
				//std::cout<<r+(lamda*v)<<"\n"<<lamda<<std::endl;
			}
			
			/*
			std::cout<<"lmin,lmax= "<<max_minus<<" "<<min_plus<<std::endl;
			std::cout<<"r= "<<r<<std::endl;
			std::cout<<"v= "<<v<<std::endl;
			std::cout<<"p1= "<<r+(min_plus*v)<<std::endl;
			std::cout<<"p2= "<<r+(max_minus*v)<<std::endl;
			
			std::pair<Point,Point> preturn(r,r);
			preturn.first[rand_coord] += min_plus;
			preturn.second[rand_coord] += max_minus;
			return preturn;	
			*/	
			
			return std::pair<NT,NT> (min_plus,max_minus);
		}
		
		std::pair<NT,NT> line_intersect_coord(Point &r,
																					Point &r_prev,
														              int rand_coord,
														              int rand_coord_prev,
														              std::vector<NT> &lamdas,
														              bool init){
			//std::cout<<"line-polytope-intersection"<<std::endl;
		  K lamda=0;
		  std::vector<NT>::iterator lamdait = lamdas.begin();
		  
		  K min_plus=0, max_minus=0;
		  bool min_plus_not_set=true;
		  bool max_minus_not_set=true;
		  
		  if(init){ //first time compute the innerprod cit*rit	
		    for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
					typename stdCoeffs::iterator cit;
					Point::Cartesian_const_iterator rit;
					rit=r.cartesian_begin(); 
					cit=ait->begin();
					K sum_nom=(*cit);
					++cit;
					K sum_denom= *(cit+rand_coord); 
					//std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
					//         v.cartesian_begin()-v.cartesian_end()<<std::endl;
					for( ; cit < ait->end() ; ++cit, ++rit){
						sum_nom -= *cit * (*rit);
					}
					lamdas[ait-_A.begin()] = sum_nom;
					if(sum_denom==K(0)){
				    //std::cout<<"div0"<<std::endl;
				    ;
				  }
				  else{
				    lamda = sum_nom*(1/sum_denom);
				    
					  if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
					  if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
					  if(lamda<min_plus && lamda>0) min_plus=lamda;
					  if(lamda>max_minus && lamda<0) max_minus=lamda;
					}
				}		
			} else {//only a few opers no innerprod
				for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
					typename stdCoeffs::iterator cit;
			    cit=ait->begin();		
					++cit;
					
					NT c_rand_coord = *(cit+rand_coord);
					NT c_rand_coord_prev = *(cit+rand_coord_prev);
					
					*lamdait = *lamdait 
						       + c_rand_coord_prev * (r_prev[rand_coord_prev] - r[rand_coord_prev]);
					
					if(c_rand_coord==K(0)){
						//std::cout<<"div0"<<std::endl;
				    ;
				  } else {					
						lamda = (*lamdait) / c_rand_coord;
						
					  if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
					  if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
					  if(lamda<min_plus && lamda>0) min_plus=lamda;
					  if(lamda>max_minus && lamda<0) max_minus=lamda;
					}
					++lamdait;
				}	  
			}	
			
			return std::pair<NT,NT> (min_plus,max_minus);
		}
		
		
	  
	private:
	  int            _d; //dimension
	  stdMatrix      _A;
};


// define different kind of polytopes
typedef std::vector<Hyperplane>           H_polytope;
typedef H_polytope 								        Polytope;
typedef std::vector<Point>						    V_polytope;

typedef std::pair<V_polytope,V_polytope> 	    MinkSumPolytope;
typedef std::pair<MinkSumPolytope,bool> 	    MinkSumPolytopeDual;


	



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

#endif //POLYTOPES_H
