// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MISC_H
#define MISC_H

//function to print rounding to double coordinates 
template <class T>
void round_print(T p) { 
    std::cout<<"test version.."<<std::endl;
   // for(typename T::Cartesian_const_iterator cit=p.cartesian_begin();
   //     cit!=p.cartesian_end(); ++cit)
    //    std::cout<<CGAL::to_double(*cit)<<" ";
   // std::cout<<std::endl;
}

/*
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
*/

// polymake file to compute exact volume
template <class T>
void print_polymake_volfile(T &P,
                            std::ostream& os){
    std::cout<<"test version.."<<std::endl;
}

/*
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
}*/

// polymake file to compute exact volume
template <class T>
void print_polymake_volfile2(T &P,
                             std::ostream& os){
    std::cout<<"test version.."<<std::endl;
}
/*
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
}*/
template <typename NT>
void read_pointset(std::istream &is,
                  std::vector<std::vector<NT> > &Input){
    
    std::string point;

    while(!std::getline(is, point, '\n').eof()) {
        //std::cout<<point<<std::endl;
        if(!std::isdigit(point[0]) && point[0]!='-' && point[0]!=' ' && point[0]!='.')
            continue;
        //std::cout<<std::endl;
        //std::size_t found = point.find_first_of(" ");
        std::size_t found =0;
        std::size_t found2=0;

        //ignore empty spaces on start of line
        found = point.find_first_not_of(" \t",found);

        std::vector<NT> input;
        while (found2!=std::string::npos || point[found]=='-')
        {
            //std::cout<<"*"<<(point[found]!='-')<<"*"<<std::endl;
            if(!std::isdigit(point[found]) && point[found]!='-')
                break;
            found2 = point.find_first_not_of("0123456789-.",found);

            //std::cout<<point.substr(found,found2-found)<<" ";
            NT num = atof(point.substr(found,found2-found).c_str());
            found=point.find_first_not_of(" \t",found2);
            //std::cout<<"found"<<point[found]<<std::endl;
            if(point[found]=='/'){
                found = found + 1;
                found2=point.find_first_not_of("0123456789-./",found);
                //std::cout<<"lala="<<point.substr(found,found2-found)<<std::endl;
                num = num / atof(point.substr(found,found2-found).c_str());
                found=point.find_first_not_of(" \t",found2);
            }
            input.push_back(num);
            //std::cout<< num<<std::endl;
        }
        Input.push_back(input);
        //std::cout<<std::endl;
    }
}

#endif //MISC_H
