#ifndef POLYTOPES_H
#define POLYTOPES_H

//#define CGAL_QP_NO_ASSERTIONS

//this is for LP-solver
#include <iostream>
//#include <CGAL/basic.h>
//#include <CGAL/QP_models.h>
//#include <CGAL/QP_functions.h>
// choose exact integral type
//#ifdef CGAL_USE_GMP
//#include <CGAL/Gmpzf.h>
//typedef CGAL::Gmpzf ET;
//#endif
//typedef double ET;
//#else
//#include <CGAL/MP_Float.h>
//typedef CGAL::MP_Float ET;
//#endif
//#include <boost/random/shuffle_order.hpp>
//#include "../rref.h"
//#include "LPsolve/solve_lp.h"
//EXPERIMENTAL
//to implement boundary oracles using NN queries  
//#include <flann/flann.hpp>



// my H-polytope class
template <typename K>
class stdHPolytope{
private:
    typedef std::vector<K>        stdCoeffs;
    typedef std::vector<stdCoeffs>  stdMatrix;
    Eigen::MatrixXd A;
    int            _d; //dimension
    stdMatrix      _A; //inequalities
    //EXPERIMENTAL
    //Flann_trees    flann_trees; //the (functional) duals of A lifted to answer NN queries
    //defined for every d coordinate
    //EXPERIMENTAL
    //typedef std::vector<flann::Index<flann::L2<double> > >  Flann_trees;

public:
    stdHPolytope() {}

    // constructor: cube(d)
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

    int dimension(){
        return _d;
    }

    int num_of_hyperplanes(){
        return _A.size();
    }

    K get_coeff(int i, int j){
        return _A[i][j];
    }

    void put_coeff(int i, int j, K value){
        _A[i][j] = value;
    }

	stdMatrix get_matrix(){
		return _A;
	}

    // default initialize: cube(d)
    int init(int d){
        A.resize(d,d);
        _d=d;
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
        return 0;
    }

    int init(stdMatrix Pin){
        _d = Pin[0][1]-1;
        
        typename stdMatrix::iterator pit=Pin.begin();
        ++pit;
        for( ; pit<Pin.end(); ++pit){
            _A.push_back(*pit);
        }
        
        //define eigen matrix
        A.resize(_A.size(),_d);
        for(int i=0; i<_A.size(); i++){
            for(int j=1; j<_d+1; j++){
                A(i,j-1)=_A[i][j];
            }
        }
        //double tstart = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::random_shuffle (_A.begin(), _A.end());
        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //std::shuffle (_A.begin(), _A.end(), std::default_random_engine(seed));
        //std::shuffle (_A.begin(), _A.end(), std::default_random_engine(seed));
        //boost::random::shuffle_order_engine<stdMatrix>(_A.begin(), _A.end());
        //double tstop = (double)clock()/(double)CLOCKS_PER_SEC;
        //std::cout << "Shuffle time = " << tstop - tstart << std::endl;
        return 0;
    }

    // print polytope in input format
    int print() {
        std::cout<<" "<<_A.size()<<" "<<_d+1<<" float"<<std::endl;
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit)
                std::cout<<*lit<<" ";
            std::cout<<std::endl;
        }
        return 0;
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

    // Compute the reduced row echelon form
    // used to transofm {Ax=b,x>=0} to {A'x'<=b'}
    // e.g. Birkhoff polytopes
    int rref(){
        to_reduced_row_echelon_form(_A);
        std::vector<int> zeros(_d+1,0);
        std::vector<int> ones(_d+1,0);
        std::vector<int> zerorow(_A.size(),0);
        for (int i = 0; i < _A.size(); ++i)
        {
            for (int j = 0; j < _d+1; ++j){
                if ( _A[i][j] == double(0)){
                    ++zeros[j];
                    ++zerorow[i];
                }
                if ( _A[i][j] == double(1)){
                    ++ones[j];
                }
            }
        }
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            int j =0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ){
                if(zeros[j]==_A.size()-1 && ones[j]==1)
                    (*mit).erase(lit);
                else{ //reverse sign in all but the first column
                    if(lit!=mit->end()-1) *lit = (-1)*(*lit);
                    ++lit;
                }
                ++j;
            }
        }
        //swap last and first columns
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            double temp=*(mit->begin());
            *(mit->begin())=*(mit->end()-1);
            *(mit->end()-1)=temp;
        }
        //delete zero rows
        for (typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ){
            int zero=0;
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit){
                if(*lit==double(0)) ++zero;
            }
            if(zero==(*mit).size())
                _A.erase(mit);
            else
                ++mit;
        }
        //update _d
        _d=(_A[0]).size();
        // add unit vectors
        for(int i=1;i<_d;++i){
            std::vector<double> e(_d,0);
            e[i]=1;
            _A.push_back(e);
        }
        // _d should equals the dimension
        _d=_d-1;
        return 1;
    }

    

    int is_in(Point p) {
        //std::cout << "Running is in" << std::endl;
        //exit(1);
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            typename stdCoeffs::iterator lit;
            //Point::Cartesian_const_iterator pit;
            //pit=p.cartesian_begin();
            typename std::vector<K>::iterator pit=p.iter_begin();
            lit=mit->begin();
            K sum=(*lit);
            ++lit;
            for( ; lit<mit->end() ; ++lit, ++pit){
                //std::cout << *lit << " " << *pit <<std::endl;
                sum -= *lit * (*pit);
            }

            //std::cout<<sum<<std::endl;
            if(sum<K(0))
                return mit-_A.begin();
        }
        return -1;
    }

    int chebyshev_center(Point& center, double& radius){
        Point f(_d);
        
        return 1;
        
    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by _A
    std::pair<Point,Point> line_intersect(Point r,
                                          Point v){
        //std::cout<<"line-polytope-intersection"<<std::endl;
        K lamda=0;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
            typename stdCoeffs::iterator cit;
            //Point::Cartesian_const_iterator rit;
            typename std::vector<K>::iterator rit=r.iter_begin();
            //rit=r.iter_begin();
            //Point::Cartesian_const_iterator vit;
            typename std::vector<K>::iterator vit=v.iter_begin();
            //vit=v.cartesian_begin();
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
        return std::pair<Point,Point> ((min_plus*v)+r,(max_minus*v)+r);
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
            typename std::vector<K>::iterator rit=r.iter_begin();
            //Point::Cartestypenameian_const_iterator rit;
            //rit=r.cartesian_begin();
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
        int mini, maxi;

        if(init){ //first time compute the innerprod cit*rit
            for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
                typename stdCoeffs::iterator cit;
                typename std::vector<K>::iterator rit=r.iter_begin();
                //Point::Cartesian_const_iterator rit;
                //rit=r.cartesian_begin();
                cit=ait->begin();
                K sum_nom=(*cit);
                ++cit;
                K sum_denom= *(cit+rand_coord);
                //std::cout<<ait->begin()-ait->end()<<" "<<r.cartesian_begin()-r.cartesian_end()<<" "<<
                //         std::endl;
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
                    //TEST
                    //if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;mini=ait-_A.begin();}
                    //if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;
                    //	 maxi=ait-_A.begin();}
                    //if(lamda<min_plus && lamda>0) {min_plus=lamda;mini=ait-_A.begin();}
                    //if(lamda>max_minus && lamda<0) {max_minus=lamda;maxi=ait-_A.begin();}
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
                    //TEST
                    //if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;mini=ait-_A.begin();}
                    //if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;
                    //	 maxi=ait-_A.begin();}
                    //if(lamda<min_plus && lamda>0) {min_plus=lamda;mini=ait-_A.begin();}
                    //if(lamda>max_minus && lamda<0) {max_minus=lamda;maxi=ait-_A.begin();}
                }
                ++lamdait;
            }
        }
        //std::cout<<"Oresult: "<<mini<<" "<<maxi<<std::endl;
        return std::pair<NT,NT> (min_plus,max_minus);
    }
    
    int linear_transformIt(Eigen::MatrixXd Tinv){
        A=A*Tinv;
        //set _A = A or replace stdMatrix with Eigen::MatrixXd ??????
        return 1;
    }

    //void rotate(){
    //std::cout<<_A<<std::endl;
    //  exit(1);
    //}


};

#endif
