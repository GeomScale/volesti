// VolEsti (volume computation and sampling library)

// Copyright (c) 2018-2020 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <Eigen/Eigen>

/// This class manipulates a point parameterized by a number type e.g. double
/// \tparam K Numerical Type
template <typename K>
class point
{
private:
    unsigned int d;

    Eigen::Matrix<typename K::FT, Eigen::Dynamic,1> coeffs;
    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    typedef Eigen::Matrix<typename K::FT, Eigen::Dynamic,1> Coeff;
    typedef typename K::FT 	FT;

    point() {}

    point(const unsigned int dim)
    {
        d = dim;
        coeffs.setZero(d);
    }

    point(const unsigned int dim, iter begin, iter endit)
    {
        d = dim;
        coeffs.resize(d);
        int i = 0;

        for (iter it=begin ; it != endit ; it++)
            coeffs(i++) = *it;
    }

    point(const Coeff& coeffs)
    {
            d = coeffs.rows();
            this->coeffs = coeffs;
    }

    point(const unsigned int dim, std::vector<typename K::FT> cofs)
    {
        d = dim;
        coeffs.resize(d);
        iter it = cofs.begin();
        int i=0;

        for (; it != cofs.end(); it++,i++)
            coeffs(i) = *it;

    }

    void add(const Coeff& coeffs)
    {
        this->coeffs += coeffs;
    }

    const Coeff& getCoefficients() const
    {
        return coeffs;
    }

    int dimension() const
    {
        return d;
    }

    void set_dimension(const unsigned int dim)
    {
        d = dim;
        coeffs.setZero(d);
    }

    void set_coord(const unsigned int i, FT coord)
    {
        coeffs(i) = coord;
    }

    void set_coeffs (const Coeff& coeffs2) {
        d = coeffs2.rows();
        coeffs = coeffs2;
    }

    void set_to_origin() {
        coeffs.setZero(d);
    }

    FT operator[] (const unsigned int i) const
    {
        return coeffs(i);
    }

    FT* pointerToData()
    {
        return coeffs.data();
    }

    FT sum() const {
        return coeffs.sum();
    }

    void operator+= (const point& p)
    {
        coeffs += p.getCoefficients();
    }

    void operator+= (const Coeff& coeffs)
    {
        this->coeffs += coeffs;
    }

    void operator-= (const point& p)
    {
        coeffs -= p.getCoefficients();
    }

    void operator-= (const Coeff& coeffs)
    {
        this->coeffs -= coeffs;
    }

    void operator= (const Coeff& coeffs)
    {
        this->coeffs = coeffs;
        d = coeffs.rows();
    }

    //TODO: avoid point construction in operators +,-,*
    point operator+ (const point& p) const
    {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs + p.getCoefficients();
        return temp;
    }

    point operator- (const point& p) const
    {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs - p.getCoefficients();
        return temp;
    }

    point operator* (const FT k) const
    {
        point temp;
        temp.d = d;
        temp.coeffs = coeffs * k;
        return temp;
    }

    void operator*= (const FT k)
    {
        coeffs *= k;
    }

    void operator/= (const FT k)
    {
        coeffs /= k;
    }

    bool operator== (point& p) const
    {
        int i=0;
        const Coeff & coeffs = p.getCoefficients();

        /* degree of approximation in
        "The art of computer programming" (vol II), p. 234, Donald. E. Knuth. */
        FT e = 0.00000000001;
        for (i=0 ; i<d ; i++) {
            if (std::abs(this->coeffs(i) - coeffs(i)) > e *std::abs(this->coeffs(i)) ||
                    std::abs(this->coeffs(i) - coeffs(i)) > e *std::abs(coeffs(i)))
                return false;
        }

        return true;
    }

    FT distance(point const & p) {
        return (this->coeffs - p.coeffs).norm();
    }

    FT dot(const point& p) const
    {
        return coeffs.dot(p.getCoefficients());
    }

    FT dot(const Coeff& coeffs) const
    {
        return this->coeffs.dot(coeffs);
    }

    FT squared_length() const {
        FT lsq = length();
        return lsq * lsq;
    }

    FT length() const {
        return coeffs.norm();
    }

    void print() const
    {
        for(unsigned int i=0; i<d; i++){
            std::cout<<coeffs(i)<<" ";
        }
        std::cout<<"\n";
    }

    static point all_ones(int dim) {
      point p(dim);
      for (int i = 0; i < dim; i++) p.set_coord(i, 1.0);
      return p;
    }

};

template<typename K>
point<K> operator* (const typename K::FT& k, point<K> const& p)
{
    return p * k;
}

#endif
