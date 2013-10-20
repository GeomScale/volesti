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

// ball type
struct Ball{
	public:
	  Ball(Point c, NT R) : _c(c),	 _R(R) {}
	   
	  Point center(){
			return _c;
		}
		NT squared_radius(){
			return _R;
		}
		NT radius(){
			return std::sqrt(_R);
		}
		int is_in(Point p){
			if ((p-_c).squared_length() <= _R)
			  return -1;
			else return 0;
		}
		std::pair<Point,Point> line_intersect(Point r, 
                                          Vector v){
			Point::Cartesian_const_iterator rit;
			rit=r.cartesian_begin(); 
			Point::Cartesian_const_iterator vit;
			vit=v.cartesian_begin();
			Point::Cartesian_const_iterator cit;
			cit=_c.cartesian_begin();
			Vector rc = r - _c;
			Vector::Cartesian_const_iterator rcit;
			rcit=rc.cartesian_begin(); 
			NT vrc(0);
			NT v2(0);
			NT rc2(0);
			for( ; cit < _c.cartesian_end() ; ++rcit, ++cit, ++rit, ++vit){
				vrc += *vit * (*rcit);
				v2 += *vit * (*vit);
				rc2 += *rcit * (*rcit);
			}
			
			//std::cout<<"r = "<<r<<std::endl;	
			//std::cout<<"v = "<<v<<std::endl;	
			//std::cout<<"_c = "<<_c<<std::endl;	
		
			//std::cout<<"vrc = "<<vrc<<std::endl;	
			//std::cout<<"v2 = "<<v2<<std::endl;	
			//std::cout<<"rc2 = "<<rc2<<std::endl;	
			//std::cout<<"disc"<<std::pow(vrc,2) - v2 * (rc2 - _R))<<std::endl;
			
			NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - _R));
			NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);									
			NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);	
			//std::cout<<"l1 = "<<lamda1<<std::endl;	
			//std::cout<<"l2 = "<<lamda2<<std::endl;	
			//std::cout<<"intersect points"<<r+(lamda1*v)<<r+(lamda2*v)<<std::endl;
			return std::pair<Point,Point> (r+(lamda1*v),r+(lamda2*v));						
		}
		
		std::pair<NT,NT> line_intersect_coord(Point r, 
                                          int rand_coord){
			
			//std::cout<<"ball intersection"<<std::endl;
			//Point::Cartesian_const_iterator rit;
			//rit=r.cartesian_begin(); 
			//Point::Cartesian_const_iterator vit;
			//vit=v.cartesian_begin();
			//Point::Cartesian_const_iterator cit;
			//cit=_c.cartesian_begin();
			//std::cout<<r<<"\n"<<rand_coord<<std::endl;
			
			Vector rc = r - _c;
			//std::cout<<"ok1"<<std::endl;
			Vector::Cartesian_const_iterator rcit;
			rcit = rc.cartesian_begin(); 
			NT vrc = *(rcit + rand_coord);
			//std::cout<<"ok2"<<std::endl;
			
			NT v2 = NT(1);
			NT rc2(0);
			for( ; rcit < rc.cartesian_end() ; ++rcit){
				//vrc += *vit * (*rcit);
				//v2 += *vit * (*vit);
				rc2 += *rcit * (*rcit);
			}
			//std::cout<<"phase I computed"<<std::endl;
			//std::cout<<"r = "<<r<<std::endl;	
			//std::cout<<"v = "<<v<<std::endl;	
			//std::cout<<"_c = "<<_c<<std::endl;	
		
			//std::cout<<"vrc = "<<vrc<<std::endl;	
			//std::cout<<"v2 = "<<v2<<std::endl;	
			//std::cout<<"rc2 = "<<rc2<<std::endl;	
			//std::cout<<"disc"<<std::pow(vrc,2) - v2 * (rc2 - _R))<<std::endl;
			
			NT disc_sqrt = std::sqrt(std::pow(vrc,2) - v2 * (rc2 - _R));
			NT lamda1((NT(-1)*vrc + disc_sqrt)/v2);									
			NT lamda2((NT(-1)*vrc - disc_sqrt)/v2);	
			//std::cout<<"l1 = "<<lamda1<<std::endl;	
			//std::cout<<"l2 = "<<lamda2<<std::endl;	
			//std::cout<<"intersect points"<<r+(lamda1*v)<<r+(lamda2*v)<<std::endl;
			
			//std::cout<<"phase II computed"<<std::endl;
			
			return std::pair<NT,NT> (lamda1,lamda2);
			
			/*
			//TODO: only change one coordinate of *r* avoid addition + construction
			std::vector<NT> v(r.dimension(),NT(0));
			v[rand_coord]=lamda1;
			Vector b1(r.dimension(),v.begin(),v.end());
			v[rand_coord]=lamda2;
			Vector b2(r.dimension(),v.begin(),v.end());
			return std::pair<Point,Point> (r+b1,r+b2);
			*/ 						
		}
	
	private:
	  Point  _c; //center
	  NT     _R; //SQUARED radius !!!
};


template <class T>
class BallIntersectPolytope {
  public:
    BallIntersectPolytope(T &P, Ball &B) : _P(P), _B(B) {};
    
    T first() { return _P; }
    Ball second() { return _B; }
    
    int is_in(Point p){
			if(_B.is_in(p)==-1)
			  return _P.is_in(p);
			return 0;
		}
		
		std::pair<Point,Point> line_intersect(Point r, 
                                          Vector v){
			
			std::pair<Point,Point> polypair = _P.line_intersect(r,v);
			std::pair<Point,Point> returnpair;
			std::pair<Point,Point> ballpair;
			bool ballinter=false;
			
			//check the first intersection point if it is inside ball
			if(_B.is_in(polypair.first)){
				//std::cout<<"inside ball 1, radius:"<<_B.radius()<<std::endl;
			  //std::cout<<polypair.first<<std::endl;
			  returnpair.first = polypair.first;
			}else{  	
				//std::cout<<"outside ball 1, radius:"<<_B.radius()<<std::endl;
				//std::cout<<polypair.first<<std::endl;
			  ballinter=true;
				//compute the intersection with ball
				ballpair = _B.line_intersect(r,v);
				returnpair.first = ballpair.first;
				//std::cout<<returnpair.first<<std::endl;
			}
			//check the second intersection point
			if(_B.is_in(polypair.second)){
			  //std::cout<<"inside ball 2, radius:"<<_B.radius()<<std::endl;
			  //std::cout<<polypair.second<<std::endl;
			  returnpair.second = polypair.second;
			}else{ 	
				//std::cout<<"outside ball 2, radius:"<<_B.radius()<<std::endl;
				//std::cout<<polypair.second<<std::endl;
			  if(ballinter) //if the intersection with ball is already computed
					returnpair.second = ballpair.second;
				else returnpair.second = (_B.line_intersect(r,v)).second;
				//std::cout<<returnpair.second<<std::endl;
			}
			return returnpair;
		}
		
		std::pair<NT,NT> line_intersect_coord(Point &r,
																					Point &r_prev,
														              int rand_coord,
														              int rand_coord_prev,
														              std::vector<NT> &lamdas,
														              bool init
														              ){
			
			std::pair<NT,NT> polypair = _P.line_intersect_coord(r,r_prev,rand_coord,rand_coord_prev,lamdas,init);
			//std::pair<NT,NT> polypair = _P.line_intersect_coord(r,r_prev,rand_coord,rand_coord_prev,lamdas,init);
			//std::cout<<"P computed"<<std::endl;
			std::pair<NT,NT> ballpair = _B.line_intersect_coord(r,rand_coord);
			//std::cout<<"P,B computed"<<std::endl;
			return std::pair<NT,NT> (std::min(polypair.first,ballpair.first),
			                         std::max(polypair.second,ballpair.second));
		}
    
  private:  
    T    _P;
    Ball _B;
};
