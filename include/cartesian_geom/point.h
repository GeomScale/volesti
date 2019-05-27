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
    unsigned int d;

    typedef Eigen::Matrix<typename K::FT, Eigen::Dynamic,1> Coeff;

    Coeff coeffs;
    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    typedef typename K::FT 	FT;

    point() {}

    point(const unsigned int dim) {
        d = dim;
        coeffs.setZero(d);
    }

    point(Coeff coeffs) {
        d = coeffs.rows();
        this->coeffs = coeffs;
    }

    point(const unsigned int dim, iter begin, iter endit) {
        d = dim;
        coeffs.resize(d);
        int i = 0;

        for (iter it=begin ; it != endit ; it++)
            coeffs(i++) = *it;
    }

    const Coeff& getCoefficients() const {
        return coeffs;
    }

    int dimension() {
        return d;
    }

    void set_dimension(const unsigned int dim) {
        d = dim;
    }

    void set_coord(const unsigned int i, FT coord) {
        coeffs(i) = coord;
    }

    FT operator[] (const unsigned int i) {
        return coeffs(i);
    }

    point operator+ (point& p) {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs + p.getCoefficients();
        return temp;
    }

    void add(Coeff coeffs) {
        this->coeffs = coeffs + this->coeffs;
    }

    point operator- (point& p) {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs - p.getCoefficients();
        return temp;
    }

    point operator* (const FT& k) {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs * k;
        return temp;
    }

    point operator/ (const FT& k) {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs / k;
        return temp;
    }

    bool operator== (point& p) {//TODO check dim?
        int i=0;

        for (auto x : p.getCoefficients()) {
            if (x != coeffs(i++)) return false;
        }

        return true;
    }


    FT dot(point& p){
        return coeffs.dot(p.getCoefficients());
    }


    FT squared_length() {

        FT lsq = FT(0.0);

        for (int i=0; i<d ; i++){
            lsq += coeffs(i) * coeffs(i);
        }
        return lsq;
    }

    void print(){
        for(unsigned int i=0; i<d; i++){
#ifdef VOLESTI_DEBUG
            std::cout<<coeffs(i)<<" ";
#endif
        }
#ifdef VOLESTI_DEBUG
        std::cout<<"\n";
#endif
    }


//    iter iter_begin() {
//        return coeffs.begin();
//    }
//
//    iter iter_end() {
//        return coeffs.end();
//    }

};

template<class K>
point<K> operator* (const typename K::FT& k, point<K>& p) {
    return p * k;
}

#endif
