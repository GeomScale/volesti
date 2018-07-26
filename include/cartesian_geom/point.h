// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POINT_H
#define POINT_H

#include <iostream>

template <class K>
class point
{
private:
    int d;
    typedef std::vector<typename K::FT> Coeff;
    Coeff coeffs;
    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    typedef typename K::FT 	FT;

    point() {}
    
    point(const int dim){
        d=dim;
        coeffs=Coeff(d,0);
    }
    
    point(const int dim, iter begin, iter end){
        d=dim;
        coeffs=Coeff(begin,end);
    }
    
    int dimension(){
        return d;
    }
    
    void set_dimension(const int dim){
        d=dim;
    }
    
    void set_coord(const int i, FT coord){
        coeffs[i]=coord;
    }
    
    FT operator[] (const int i){
        return coeffs[i];
    }
    
    point operator+ (point& p){
        point temp(p.dimension());
        typename Coeff::iterator tmit=temp.iter_begin();
        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit, ++tmit){
			(*tmit)=(*pit)+(*mit);
		}
        //for (int i=0; i<d; i++){
            //temp.coeffs[i]=p[i]+coeffs[i];
        //}
        return temp;
    }
    
    point operator- (point& p){
        point temp(p.dimension());
        typename Coeff::iterator tmit=temp.iter_begin();
        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit, ++tmit){
			(*tmit)=(*mit)-(*pit);
		}
        //for (int i=0; i<d; i++){
            //temp.coeffs[i]=coeffs[i]-p[i];
        //}
        
        return temp;
    }
    
    point operator* (const FT& k){
        //point temp(d,iter_begin(),iter_end());
        point temp(d);
        //point temp=*this;
        //typename Coeff::iterator tmit=temp.iter_begin();
        //for( ; tmit<temp.iter_end(); ++tmit){
		//	(*tmit)*=k;
		//}
        for (int i=0; i<d; i++){
            temp.coeffs[i]=coeffs[i]*k;
        }
        
        return temp;
    }


    FT dot(point& p){
        FT res=FT(0);

        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit){
            res+=(*mit)*(*pit);
        }
        return res;
    }
    
    
    FT squared_length(){
        
        FT lsq=FT(0);
        
        for (int i=0; i<d; i++){
            lsq+=coeffs[i]*coeffs[i];
        }    
        return lsq;
    }

    void print(){
        std::cout<<"\n";
        for(int i=0; i<d; i++){
            std::cout<<coeffs[i]<<" ";
        }
        std::cout<<"\n";
    }
    
    
    iter iter_begin(){
        
        return coeffs.begin();
        
    }
    
    iter iter_end(){
        
        return coeffs.end();
        
    }
    
    
};

template<class K>
point<K> operator* (const typename K::FT& k, point<K>& p){
        
        return p*k;
        
}

#endif
