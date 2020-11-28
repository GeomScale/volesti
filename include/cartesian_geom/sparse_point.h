// VolEsti (volume computation and sampling library)

// Copyright (c) 2018-2020 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SPARSE_POINT_H
#define SPARSE_POINT_H

#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

template <typename K>
class sparse_point
{
private:
    unsigned int d;

    typedef typename std::vector<typename K::FT>::iterator iter;
public:
    Eigen::SparseMatrix<typename K::FT> coeffs;
    typedef Eigen::SparseMatrix<typename K::FT> Coeff;
    typedef typename K::FT 	FT;

    sparse_point() {}

    sparse_point(const unsigned int dim)
    {
        d = dim;
        Coeff coeffs(d, 1);
    }

    sparse_point(const Coeff& coeffs)
    {
            d = coeffs.rows();
            this->coeffs = coeffs;
    }

    sparse_point(const unsigned int dim, std::vector<typename K::FT> cofs)
    {
        d = dim;
        coeffs = Coeff(d);
        coeffs.resize(d);
        iter it = cofs.begin();
        int i=0;

        for (; it != cofs.end(); it++,i++)
            coeffs.coeffRef(i, 0) = *it;

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
    }

    void set_coord(const unsigned int i, FT coord)
    {
        coeffs.coeffRef(i, 0) = coord;
    }

    void set_coeffs (const Coeff& coeffs2) {
        d = coeffs2.rows();
        coeffs = coeffs2;
    }

    void set_to_origin() {
        coeffs.setZero();
    }

    FT operator[] (const unsigned int i) const
    {
        return coeffs(i);
    }

    FT* sparse_pointerToData()
    {
        return coeffs.data();
    }

    FT sum() const {
        return coeffs.sum();
    }

    void operator+= (const sparse_point& p)
    {
        coeffs += p.getCoefficients();
    }

    void operator+= (const Coeff& coeffs)
    {
        this->coeffs = coeffs + this->coeffs;
    }

    void operator= (const Coeff& coeffs)
    {
        this->coeffs = coeffs;
        d = coeffs.rows();
    }

    //TODO: avoid sparse_point construction in operators +,-,*
    sparse_point operator+ (const sparse_point& p) const
    {
        sparse_point temp;
        temp.d = d;
        temp.coeffs = coeffs + p.getCoefficients();
        return temp;
    }

    sparse_point operator- (const sparse_point& p) const
    {
        sparse_point temp;
        temp.d = d;
        temp.coeffs = coeffs - p.getCoefficients();
        return temp;
    }

    sparse_point operator* (const FT k) const
    {
        sparse_point temp;
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

    bool operator== (sparse_point& p) const
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

    FT distance(sparse_point const & p) {
        return (this->coeffs - p.coeffs).norm();
    }

    FT dot(const sparse_point& p) const
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
        for(unsigned int i =  0; i<d; i++){
            std::cout<< coeffs.coeff(i, 0)<<" ";
        }
        std::cout<<"\n";
    }

    static sparse_point all_ones(int dim) {
      sparse_point p(dim);
      for (int i = 0; i < dim; i++) p.coeffs.coeffRef(i, 0) = FT(1.0);
      return p;
    }

};

template<typename K>
sparse_point<K> operator* (const typename K::FT& k, sparse_point<K> const& p)
{
    return p * k;
}

#endif
