// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <cmath>
#include <vector>

template <class K>
class point
{
private:
    unsigned int d;
    typedef typename std::vector<typename K::FT> Coeff;
    Coeff coeffs;
    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    typedef typename K::FT 	FT;

    point() {}
    
    point(const unsigned int dim) {
        d = dim;
        coeffs = Coeff(d,0);
    }
    
    point(const unsigned int dim, iter begin, iter endit) {
        d = dim;
        coeffs = Coeff(begin,endit);
    }

    point(const unsigned int dim, FT* data) {
        d = dim;
        coeffs.reserve(d);
        for (uint i=0; i<dim; i++) {
            coeffs[i] = data[i];
        }
    }

    FT* data() {
        return coeffs.data();
    }

    Coeff& get_coeffs() {
        return coeffs;
    }
    
    int dimension() {
        return d;
    }
    
    void set_dimension(const unsigned int dim) {
        d = dim;
        coeffs.reserve(d);
    }
    
    void set_coord(const unsigned int i, FT coord) {
        coeffs[i] = coord;
    }
    
    FT operator[] (const unsigned int i) {
        return coeffs[i];
    }
    
    point operator+ (point& p) {
        point temp(p.dimension());

        typename Coeff::iterator tmit = temp.iter_begin();
        typename Coeff::iterator pit = p.iter_begin();
        typename Coeff::iterator mit = coeffs.begin();

        for (; pit < p.iter_end(); ++pit, ++mit, ++tmit) {
            (*tmit) = (*pit) + (*mit);
        }
        return temp;
    }
    
    point operator- (point& p) {
        point temp(p.dimension());

        typename Coeff::iterator tmit = temp.iter_begin();
        typename Coeff::iterator pit = p.iter_begin();
        typename Coeff::iterator mit = coeffs.begin();

        for (; pit < p.iter_end(); ++pit, ++mit, ++tmit) {
            (*tmit) = (*mit) - (*pit);
        }
        return temp;
    }

    double squared_distance(point& p) {

        double dist = 0.0;
        for (auto mit = coeffs.begin(), pit=p.iter_begin(); pit < p.iter_end(); ++pit, ++mit) {
            double tmp = (*mit)-(*pit);
            dist += tmp*tmp;
        }
        return dist;
    }

    point operator* (const FT& k) {
        point temp(d, iter_begin(), iter_end());

        typename Coeff::iterator tmit = temp.iter_begin();

        for (; tmit < temp.iter_end(); ++tmit) {
            (*tmit) *= k;
        }
        return temp;
    }


    FT dot(point& p){
        FT res=FT(0.0);

        typename Coeff::iterator pit=p.iter_begin();
        typename Coeff::iterator mit=coeffs.begin();
        for( ; pit<p.iter_end(); ++pit, ++mit){
            res+=(*mit)*(*pit);
        }
        return res;
    }
    
    
    FT squared_length() {

        FT lsq = FT(0.0);

        typename Coeff::iterator mit=coeffs.begin();
        for ( ; mit != coeffs.end(); mit++){
            lsq += (*mit)*(*mit);
        }
        return lsq;
    }

    void print(){
        for(unsigned int i=0; i<d; i++){
            #ifdef VOLESTI_DEBUG
            std::cout<<coeffs[i]<<" ";
            #endif
        }
        #ifdef VOLESTI_DEBUG
        std::cout<<"\n";
        #endif
    }

    void normalize() {
        FT lsq=std::sqrt(squared_length());
        for (auto mit=coeffs.begin(); mit!=coeffs.end(); mit++){
            (*mit) = (*mit)/lsq;
        }
    }
    
    iter iter_begin() {
        return coeffs.begin();
    }
    
    iter iter_end() {
        return coeffs.end();
    }
    
    
};

template<class K>
point<K> operator* (const typename K::FT& k, point<K>& p) {
    return p * k;
}

#endif
