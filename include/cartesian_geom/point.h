// VolEsti (volume computation and sampling library)

// Copyright (c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POINT_H
#define POINT_H

#include <iostream>
//#include "Eigen"


template <class K>
class point
{
private:
    unsigned int d;

    typedef typename K::VT Coeff;
    typedef typename K::MT MT;
    //typedef Eigen::Matrix<typename K::FT,Eigen::Dynamic,Eigen::Dynamic> MT;

    Coeff coeffs;
    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    typedef typename K::FT 	FT;

    point() {}

    point(const unsigned int& dim) {
        d = dim;
        coeffs.setZero(d);
    }

    point(const Coeff& coeffs) {
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

    int dimension() const {
        return d;
    }

    void set_dimension(const unsigned int dim) {
        d = dim;
    }

    void set_coord(const unsigned int i, FT coord) {
        coeffs(i) = coord;
    }

    void set(const Coeff& vector) {
        this->coeffs.noalias() = vector;
    }

    Coeff get_coefficients() {
        return coeffs;
    }

    FT operator[] (const unsigned int& i) const {
        return coeffs(i);
    }

    point operator+ (const point& p) const {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs + p.getCoefficients();
        return temp;
    }

    void add(const Coeff& coeffs) {
        this->coeffs = coeffs + this->coeffs;
    }

    void operator+= (const point& p)  {
        this->coeffs += p.coeffs;
    }

    void operator*= (const FT& k)  {
        this->coeffs *= k;
    }

    void operator/= (const FT& k)  {
        this->coeffs /= k;
    }

    void operator-= (const point& p)  {
        this->coeffs -= p.coeffs;
    }

    point operator- (const point& p) const{
        point temp;
        temp.d = d;
        temp.coeffs = coeffs - p.getCoefficients();
        return temp;
    }

    point operator* (const FT& k) const {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs * k;
        return temp;
    }

    point operator/ (const FT& k) const{
        point temp;
        temp.d = d;
        temp.coeffs = coeffs / k;
        return temp;
    }

    bool operator== (const point& p) const{//TODO check dim?
        int i=0;

        /* degree of approximation in
        "The art of computer programming" (vol II), p. 234, Donald. E. Knuth. */
        FT e = 0.00000000001;

        for (i=0 ; i<d ; i++) {
            if (std::abs(coeffs(i) - p[i]) > e * std::abs(p[i]) ||
                std::abs(coeffs(i) - p[i]) > e * std::abs(coeffs(i))){
                return false;
            }
        }

        return true;
    }


    FT dot(const point& p) const{
        return coeffs.dot(p.getCoefficients());
    }

    point matrix_left_product(const MT& matrix) const {
        return point(matrix* coeffs);
    }

    double normalize() {
        double norm = coeffs.norm();
        this->coeffs.normalize();
        return norm;
    }

    FT squared_length() {

        FT lsq = FT(0.0);

        for (int i=0; i<d ; i++){
            lsq += coeffs(i) * coeffs(i);
        }
        return lsq;
        // TODO: return std::pow(coeffs.norm(), FT(2.0));
    }

    void print(){
        for(unsigned int i=0; i<d; i++){
//#ifdef VOLESTI_DEBUG
            std::cout<<coeffs(i)<<" ";
//#endif
        }
//#ifdef VOLESTI_DEBUG
        std::cout<<"\n";
//#endif
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
point<K> operator* (const typename K::FT& k, const point<K>& p) {
    return p * k;
}

#endif
