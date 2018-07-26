// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <iostream>


// my H-polytope class
template <typename K>
class Polytope{
private:
    typedef std::vector<K>        stdCoeffs;
    typedef std::vector<stdCoeffs>  stdMatrix;
    Eigen::MatrixXd A;
    int            _d; //dimension
    stdMatrix      _A; //inequalities

public:
    typedef K                    FT;
    Polytope() {}

    // constructor: cube(d)
    Polytope(int d): _d(d) {
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
        //print();
        //std::cout<<"A eigen: \n"<<A<<std::endl;
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
        for(typename stdMatrix::iterator mit=_A.begin(); mit<_A.end(); ++mit){
            typename stdCoeffs::iterator lit;
            typename std::vector<K>::iterator pit=p.iter_begin();
            lit=mit->begin();
            K sum=(*lit);
            ++lit;
            for( ; lit<mit->end() ; ++lit, ++pit){
                sum -= *lit * (*pit);
            }
            if(sum<K(0))
                return mit-_A.begin();
        }
        return -1;
    }

    std::pair<Point,double> chebyshev_center(){

        std::pair<Point,double> res;
        res=solveLP(_A,_d);
        return res;
        
    }

    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by _A
    /*
    std::pair<Point,Point> line_intersect(Point r,
                                          Point v){
        K lamda=0;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
            typename stdCoeffs::iterator cit;
            typename std::vector<K>::iterator rit=r.iter_begin();
            typename std::vector<K>::iterator vit=v.iter_begin();
            cit=ait->begin();
            K sum_nom=(*cit);
            ++cit;
            K sum_denom=K(0);
            for( ; cit < ait->end() ; ++cit, ++rit, ++vit){
                sum_nom -= *cit * (*rit);
                sum_denom += *cit * (*vit);
            }
            if(sum_denom==K(0)){
                //std::cout<<"div0"<<std::endl;
                ;
            }
            else{
                lamda = sum_nom/sum_denom;
                if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
                if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
                if(lamda<min_plus && lamda>0) min_plus=lamda;
                if(lamda>max_minus && lamda<0) max_minus=lamda;
            }
        }
        
        return std::pair<Point,Point> ((min_plus*v)+r,(max_minus*v)+r);
    }
     */
    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by _A
    std::pair<NT,NT> line_intersect(Point r,
                                          Point v){
        K lamda=0;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
            typename stdCoeffs::iterator cit;
            typename std::vector<K>::iterator rit=r.iter_begin();
            typename std::vector<K>::iterator vit=v.iter_begin();
            cit=ait->begin();
            K sum_nom=(*cit);
            ++cit;
            K sum_denom=K(0);
            for( ; cit < ait->end() ; ++cit, ++rit, ++vit){
                sum_nom -= *cit * (*rit);
                sum_denom += *cit * (*vit);
            }
            if(sum_denom==K(0)){
                //std::cout<<"div0"<<std::endl;
                ;
            }
            else{
                lamda = sum_nom/sum_denom;
                if(min_plus_not_set && lamda>0){min_plus=lamda;min_plus_not_set=false;}
                if(max_minus_not_set && lamda<0){max_minus=lamda;max_minus_not_set=false;}
                if(lamda<min_plus && lamda>0) min_plus=lamda;
                if(lamda>max_minus && lamda<0) max_minus=lamda;
            }
        }

        //return std::pair<Point,Point> ((min_plus*v)+r,(max_minus*v)+r);
        return std::pair<NT,NT> (min_plus, max_minus);
    }


    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          int rand_coord){
        K lamda=0;
        K min_plus=0, max_minus=0;
        bool min_plus_not_set=true;
        bool max_minus_not_set=true;
        for(typename stdMatrix::iterator ait=_A.begin(); ait<_A.end(); ++ait){
            typename stdCoeffs::iterator cit;
            typename std::vector<K>::iterator rit=r.iter_begin();
            cit=ait->begin();
            K sum_nom=(*cit);
            ++cit;
            K sum_denom= *(cit+rand_coord);
            for( ; cit < ait->end() ; ++cit, ++rit){
                sum_nom -= *cit * (*rit);
            }
            if(sum_denom==K(0)){
                //std::cout<<"div0"<<sum_denom<<std::endl;
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
        return std::pair<NT,NT> (min_plus,max_minus);
    }

    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          int rand_coord,
                                          int rand_coord_prev,
                                          std::vector<NT> &lamdas,
                                          bool init){
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
                cit=ait->begin();
                K sum_nom=(*cit);
                ++cit;
                K sum_denom= *(cit+rand_coord);
                for( ; cit < ait->end() ; ++cit, ++rit){
                    sum_nom -= *cit * (*rit);
                }
                lamdas[ait-_A.begin()] = sum_nom;
                if(sum_denom==K(0)){
                    //std::cout<<"div0"<<sum_denom<<std::endl;
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
        //std::cout<<"Oresult: "<<mini<<" "<<maxi<<std::endl;
        return std::pair<NT,NT> (min_plus,max_minus);
    }
    
    int linear_transformIt(Eigen::MatrixXd Tinv){
        A=A*Tinv;
        //set _A = A or replace stdMatrix with Eigen::MatrixXd ??????
        for(int i=0; i<_A.size(); i++){
            for(int j=1; j<_d+1; j++){
                _A[i][j]=A(i,j-1);
            }
        }
        return 1;
    }



};

#endif
