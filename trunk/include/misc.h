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

#ifndef MISC_H
#define MISC_H

//function to print rounding to double coordinates 
template <class T>
void round_print(T p) { 
  for(typename T::Cartesian_const_iterator cit=p.cartesian_begin(); 
      cit!=p.cartesian_end(); ++cit)
	  std::cout<<CGAL::to_double(*cit)<<" "; 
  std::cout<<std::endl;
}

// Naive algorithm for Mink sum
typedef std::vector<V_polytope>              Vpolys;

int Minkowski_sum_naive(V_polytope &P1, V_polytope &P2, V_polytope &Msum){
	std::cout<<(!P1.empty() && !P2.empty())<<std::endl;
	if(!P1.empty() && !P2.empty()){
	  V_polytope Msum_all;
		for (V_polytope::iterator Pit1 = P1.begin(); Pit1 != P1.end(); ++Pit1){
	    for (V_polytope::iterator Pit2 = P2.begin(); Pit2 != P2.end(); ++Pit2){
	      Point p = CGAL::Origin() + 
	            (((*Pit1)-CGAL::Origin()) + ((*Pit2)-CGAL::Origin()));
	      Msum_all.push_back(p);
	      //std::cout<<p<<std::endl;
	    }
	  } 
	  //std::cout<<"---------"<<std::endl;
	  // compute the extreme points
		CGAL::Extreme_points_d<EP_Traits_d> ep(P1[0].dimension());
	  ep.insert(Msum_all.begin(),Msum_all.end());
		//std::vector<Point> extreme_points;
		ep.get_extreme_points(std::back_inserter(Msum));
		return Msum.size();
  }
  return -1;
}


// polymake file to compute exact volume
template <class T>
void print_polymake_volfile(T &P,
														std::ostream& os){
  // print the vertices of the P polytope
  os << "use Time::HiRes qw(gettimeofday tv_interval);\n";
	os << "use application 'polytope';\n";
	os << "my $p=new Polytope<Rational>;\n";
	os << "$p->POINTS=<<'.';\n";
  for (typename T::iterator vit = P.begin(); vit != P.end(); vit++){
    os << "1 ";
    for (Point::Cartesian_const_iterator cit=vit->cartesian_begin();
         cit != vit->cartesian_end();
         cit++){
      os << *cit;
      if (cit - vit->cartesian_begin() != vit->dimension()-1)
        os << " ";
    }
    //os << "|" << vit->point().index();
    os << "\n";
  }
	os << ".\n";
	os << "print ' ';\n";
	os << "print $p->N_POINTS;\n";
	os << "print ' ';\n";
	os << "print $p->N_VERTICES;\n";
	os << "print ' ';\n";
	os << "print $p->DIM;\n";
	os << "print ' ';\n";
	os << "my $t0 = [gettimeofday];\n";
	os << "my $f=$p->VOLUME;\n";
	os << "print $f;\n";
	os << "print ' ';\n";
	os << "print tv_interval($t0,[gettimeofday]);\n";
	os << "print \"\n\";\n";
	os << std::endl;
}

// polymake file to compute exact volume
template <class T>
void print_polymake_volfile2(T &P,
														std::ostream& os){
  // print the vertices of the P polytope
  os << "use Time::HiRes qw(gettimeofday tv_interval);\n";
	os << "use application 'polytope';\n";
	os << "my $p=new Polytope;\n";
	os << "$p->INEQUALITIES=<<'.';\n";
	//os << "my $p=new Polytope<Rational>;\n";
	//os << "$p->POINTS=<<'.';\n";
  for (typename T::iterator vit = P.begin(); vit != P.end(); vit++){
    Hyperplane::Coefficient_const_iterator cit_end = vit->coefficients_end();
    os << *(--cit_end)<<" ";
    //os << "0 ";
    Hyperplane::Coefficient_const_iterator cit = vit->coefficients_begin();
    //++cit;
    for (; cit != cit_end; cit++){
      //std::cout<<*cit<<" ";
      os <<(*cit)<<" ";
      if (cit - vit->coefficients_begin() != vit->dimension()-1)
        os << " ";
    }
    //os << "|" << vit->point().index();
    os << "\n";
  }
	os << ".\n";
	//$p=new Polytope<Rational>(INEQUALITIES=>$inequalities);

	os << "print ' ';\n";
	os << "print $p->N_POINTS;\n";
	os << "print ' ';\n";
	os << "print $p->N_VERTICES;\n";
	os << "print ' ';\n";
	os << "print $p->DIM;\n";
	os << "print ' ';\n";
	os << "my $t0 = [gettimeofday];\n";
	os << "my $f=$p->VOLUME;\n";
	os << "print $f;\n";
	os << "print ' ';\n";
	os << "print tv_interval($t0,[gettimeofday]);\n";
	os << "print \"\n\";\n";
	os << std::endl;
}



#endif //MISC_H
