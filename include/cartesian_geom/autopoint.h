// VolEsti (volume computation and sampling library)

// Copyright (c) 2024 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef CARTESIAN_KERNEL_AUTOPOINT_H
#define CARTESIAN_KERNEL_AUTOPOINT_H

#include <iostream>
#include <Eigen/Eigen>

/// This class manipulates a point used for automatic differentation
/// parameterized by a number type e.g. double
/// \tparam T Numerical Type
template <typename T>
class autopoint
{
    public:
    unsigned int d;
    Eigen::Matrix<typename autodiff::detail::Real<1, T>, Eigen::Dynamic,1> coeffs;
    typedef typename std::vector<typename autodiff::detail::Real<1, T>>::iterator iter;
    typedef typename autodiff::detail::Real<1, T> 	FT;
    typedef Eigen::Matrix<T, Eigen::Dynamic,1> coeff;
    typedef Eigen::Matrix<typename autodiff::detail::Real<1, T>, Eigen::Dynamic,1> Coeff;
    autopoint() {}

    autopoint(const unsigned int dim)
    {
        d = dim;
        coeffs.setZero(d);
    }

    autopoint(const unsigned int dim, iter begin, iter endit)
    {
        d = dim;
        coeffs.resize(d);
        int i = 0;

        for (iter it=begin ; it != endit ; it++)
            coeffs(i++) = *it;
    }
    FT operator()(int i,int j)
    {
        return coeffs(i,j);
    }

    autopoint(const Coeff& coeffs)
    {
            d = coeffs.rows();
            this->coeffs = coeffs;
    }

    autopoint(const coeff& coeffs)
    {
        d = coeffs.rows();
        // Coeff temp= coeffs;
        this->coeffs = (Coeff)coeffs;
    }


    autopoint(const unsigned int dim, std::vector<typename autodiff::detail::Real<1, T>> cofs)
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
    autopoint head(int n) const
    {
        return autopoint((Coeff)coeffs.head(n));
    }

        autopoint tail(int n) const
    {
        return autopoint((Coeff)coeffs.tail(n));
    }

    autopoint pow(int n) const {
        autopoint temp;
        temp.d = d;
        temp.coeffs=coeffs.array().pow(n).matrix();
        return temp;
    }


    autopoint log() const {
        autopoint temp;
        temp.d = d;
        temp.coeffs=coeffs.array().log().matrix();
        return temp;
    }

    autopoint exp() const {
        autopoint temp;
        temp.d = d;
        temp.coeffs=coeffs.array().exp().matrix();
        return temp;
    }
    void operator+= (const autopoint& p)
    {
        coeffs += p.getCoefficients();
    }

    void operator+= (const Coeff& coeffs)
    {
        this->coeffs += coeffs ;
    }

    void operator-= (const autopoint& p)
    {
        coeffs -= p.getCoefficients();
    }

    void operator-= (const Coeff& coeffs)
    {
        this->coeffs -= coeffs ;
    }


    void operator= (const Coeff& coeffs)
    {
        this->coeffs = coeffs;
        d = coeffs.rows();
    }

    // copy assignment
    autopoint& operator=(const autopoint& coeffs)
    {
    // Guard self assignment
        if (this == &coeffs)
            return *this;
        this->coeffs = coeffs;
        d = coeffs.rows();
        return *this;
    }
    autopoint transpose() const
    {
        return autopoint((Coeff)coeffs.transpose());
    }

    autopoint operator+ (const autopoint& p) const
    {
        return autopoint((Coeff)(coeffs+p.getCoefficients()));
    }
        autopoint operator- (const autopoint& p) const
    {
        return autopoint((Coeff)(coeffs-p.getCoefficients()));
    }

    autopoint operator- (const FT p) const
    {
         autopoint temp_auto;
         temp_auto.d=d;
         auto temp=coeffs.array();
         temp_auto.coeffs=(temp-p).matrix();
         return temp_auto;
    }
    autopoint operator- (T p) const
    {
         autopoint temp_auto;
         temp_auto.d=d;
         auto temp=coeffs.array();
         temp_auto.coeffs=(temp-p).matrix();
         return temp_auto;
    }

    autopoint operator* (const FT k)
    {
        return autopoint(coeffs * k);
    }

    autopoint operator* (T k) const
    {
        return autopoint((Coeff)(coeffs * ((FT)k)));
    }

    autopoint operator* (T k)
    {
        return autopoint((Coeff)(coeffs * ((FT)k)));
    }
    // matrix multiplication
    autopoint operator* (const autopoint& autopoint_)
    {
        return autopoint((Coeff)(coeffs * autopoint_.getCoefficients()));
    }
    // matrix multiplication with normal matrix
    autopoint operator* (const coeff& matrix_)
    {
        return autopoint((  autopoint(matrix_) * autopoint(coeffs) ));
    }

    void operator*= (const FT k)
    {
        coeffs *= k;
    }

    void operator*= (const T k)
    {
        FT k_=k;
        coeffs *= k_;
    }


    void operator/= (const FT k)
    {
        coeffs /= k;
    }


    FT distance(const autopoint& p) const
    {
        return (this->coeffs - p.coeffs).norm();
    }

    FT dot(const autopoint& p) const
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

    static autopoint all_ones(int dim)
    {
      autopoint p(dim);
      for (int i = 0; i < dim; i++) p.set_coord(i, 1.0);
      return p;
    }

};

template<typename K>
autopoint<K> operator* ( K k, autopoint<K> const& p)
{
    return p * k;
}

// matrix times autopoint
template<typename K>
autopoint<K> operator* (  Eigen::Matrix<K, Eigen::Dynamic,Eigen::Dynamic> matrix_, autopoint<K> const& p)
{
    Eigen::Matrix<typename autodiff::detail::Real<1, K>, Eigen::Dynamic,1> temp= matrix_*p.getCoefficients();
    return autopoint(temp);
}



#endif // CARTESIAN_KERNEL_AUTOPOINT_H
