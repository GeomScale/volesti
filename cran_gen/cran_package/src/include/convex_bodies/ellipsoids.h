// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

#include <iostream>


template <class Point, class MT, class VT>
class copula_ellipsoid{
private:
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    //typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    //typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    MT G;
    unsigned int dim;
public:

    copula_ellipsoid() {}

    copula_ellipsoid(std::vector<std::vector<NT> > Gin) {
        dim = Gin.size();
        G.resize(dim, dim);
        for (unsigned int i = 0; i < Gin.size(); i++) {
            for (unsigned int j = 0; j < Gin.size(); j++) {
                G(i,j) = Gin[i][j];
            }
        }
    }

    NT mat_mult(Point p) {
         return p.getCoefficients().transpose()*G*p.getCoefficients();
    }

};


/* developing part
// ellipsoid class
template <typename K>
class Ellipsoid{
private:
    typedef std::vector<K>        stdCoeffs;
    typedef std::vector<stdCoeffs>  stdMatrix;
    int d;
    stdMatrix C;
    K c0;
    
public:
    Ellipsoid(){}

    Ellipsoid(int dim, stdMatrix Cin, K c0in){
        d=dim;
        typename stdMatrix::iterator pit=Cin.begin();
        for( ; pit<Cin.end(); ++pit){
            C.push_back(*pit);
        }
    }
    
    int dimension(){
        return d;
    }
    
    K get_coeff(int i, int j){
        return C[i][j];
    }

    void put_coeff(int i, int j, K value){
        C[i][j] = value;
    }
    
    void print() {
        #ifdef VOLESTI_DEBUG
        std::cout<<" "<<C.size()<<" "<<d+1<<" float"<<std::endl;
        #endif
        for(typename stdMatrix::iterator mit=C.begin(); mit<C.end(); ++mit){
            for(typename stdCoeffs::iterator lit=mit->begin(); lit<mit->end() ; ++lit){
                #ifdef VOLESTI_DEBUG
                std::cout<<*lit<<" ";
                #endif
            }
            #ifdef VOLESTI_DEBUG
            std::cout<<std::endl;
            #endif
        }
    }
    
    int is_in(Point p){
        K sum=K(0);
        for(typename stdMatrix::iterator cit=C.begin(); cit<C.end(); ++cit){
            typename stdCoeffs::iterator rit;
            //Point::Cartesian_const_iterator pit;
            //pit=p.cartesian_begin();
            typename std::vector<K>::iterator pit=p.iter_begin();
            rit=cit->begin();
            //K sum=(*lit);
            //++lit;
            for( ; rit<cit->end() ; ++rit, ++pit){
                sum += (*rit) * (*pit);
            }

            //std::cout<<sum<<std::endl;
        }
        if(sum>c0){
            return 0;
        }
        return -1;
        
    }
    
    std::pair<Point,Point> line_intersect(Point p,
                                          Point v){
        K a=K(0) , b1=K(0) , b2=K(0) , c=K(0) , D;
        stdCoeffs Cu(d,K(0)), Cx(d,K(0));
        int i;
        
        for(typename stdMatrix::iterator cit=C.begin(); cit<C.end(); ++cit){
            typename stdCoeffs::iterator rit;
            typename std::vector<K>::iterator vit=v.iter_begin();
            typename std::vector<K>::iterator pit=p.iter_begin();
            rit=cit->begin();
            i=0;
            for( ; rit<cit->end() ; ++rit, ++pit, ++vit, ++i){
                Cu[i]+=(*rit) * (*vit);
                Cx[i]+=(*rit) * (*pit);
            }
        }
        for (i=0; i<d; i++){
            a+=v[i]*Cu[i];
            b1+=v[i]*Cx[i];
            b2+=p[i]*Cu[i];
            c+=p[i]*Cx[i];
        }
        b1+=b2;
        
        D=std::pow(b1,2)-4*a*c;
        return std::pair<Point,Point> ( ( ((-b1+std::sqrt(D))/(2*a))*v)+p , ( ((-b1-std::sqrt(D))/(2*a))*v)+p );
    }
    
    std::pair<K,K> line_intersect_coord(Point &p,
                                          int rand_coord){
        K a=K(0) , b=K(0) , c=K(0) , D;
        int i,j;

	    for (i=0; i<d; i++) {
            if (i == rand_coord) {
                b += 2 * C[i][i] * p[i];
                a += C[i][i];
            }
            c += C[i][i] * std::pow(p[i], 2);
            for (j = i + 1; j < d; j++) {
                if (i == rand_coord) {
                    b += 2 * C[i][j] * p[j];
                }
                if (j == rand_coord) {
                    b += 2 * C[i][j] * p[i];
                }
                c += C[i][j] * p[i] * p[j] * 2;
            }
        }
        c-=c0;
        D=std::pow(b,2)-4*a*c;
        return std::pair<K,K> ((-b+std::sqrt(D))/(2*a) , (-b-std::sqrt(D))/(2*a));

    }
    
}; */

#endif
