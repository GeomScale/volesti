// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// This file is converted from BNMin1 (https://www.mrao.cam.ac.uk/~bn204/oof/bnmin1.html) by Apostolos Chalkis

// Original copyright notice:

/**
   Bojan Nikolic <bojan@bnikolic.co.uk>
   Initial version 2010

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file ellipsoids.cxx

   Computation and use of ellipsoids releated to sets of points
*/
#ifndef KHACH_H
#define KHACH_H

#include <set>
#include <iostream>
#include <Eigen/Eigen>

//#include "khach1.h"
//#include "mcpoint1.h"
#include "mcpoint.h"
//#include "bnmin_main1.h"
//#include "bnmin_main2.h"

//#include "../bnmin_main.hxx"

//namespace Minim {

  template <class NT>
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;

  template <class NT>
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

  struct KhachiyanEllipsoid
  {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Q;
      Eigen::Matrix<double, Eigen::Dynamic, 1>              c;
  };

  template<typename Derived>
  inline bool is_nan(const Eigen::MatrixBase<Derived>& x)
  {
  return ((x.array() == x.array())).all();
  }

  template<class T>
  bool InvertMatrix(const MT<T> &input,
                    MT<T> &inverse)
  {
    inverse = input.inverse();
    return !is_nan(inverse);
    // return true;
  }


  inline void InvertLP(const MT<double> &Lambdap,
                MT<double> &LpInv)
  {
    bool res = InvertMatrix(Lambdap, LpInv);
    if (not res)
    {
     // throw an error of your choice here!
     // throw MatrixErr("Could not invert matrix: ",
     //                 Lambdap);
    }
  }

  inline void Lift(const MT<double> &A, MT<double> &Ap)
  {
    Ap.resize(A.rows()+1, A.cols());
    Ap.topLeftCorner(A.rows(), A.cols()) = A;
    Ap.row(Ap.rows()-1).setConstant(1.0); 
  }

  inline void genDiag(const VT<double> &p, MT<double> &res)
  {
    res.setZero();

    for(size_t i=0; i<p.rows(); ++i)
    {
      res(i,i)=p(i);
    }
  }

  inline void KaLambda(const MT<double> &Ap,
                const VT<double> &p,
                MT<double> &Lambdap)
  {

    MT<double> dp(p.size(), p.size());
    genDiag(p, dp);

    dp = dp * Ap.transpose();
    Lambdap.noalias() = Ap * dp;
  }

  inline double KhachiyanIter(const MT<double> &Ap, VT<double> &p)
  {
    /// Dimensionality of the problem
    const size_t d = Ap.rows()-1;

    MT<double> Lp;
    MT<double> M;
    KaLambda(Ap, p, Lp);
    MT<double> ILp(Lp.rows(), Lp.cols());
    InvertLP(Lp, ILp);
    M.noalias() = ILp * Ap;
    M = Ap.transpose() * M;

    double maxval=0;
    size_t maxi=0;
    for(size_t i=0; i<M.rows(); ++i)
    {
      if (M(i,i) > maxval)
      {
        maxval=M(i,i);
        maxi=i;
      }
    }
    const double step_size=(maxval -d - 1)/((d+1)*(maxval-1));
    VT<double> newp = p*(1-step_size);
    newp(maxi) += step_size;

    const double err= (newp-p).norm();
    p = newp;
    return err;

  }

  inline void KaInvertDual(const MT<double> &A, 
                      const VT<double> &p, 
                      MT<double> &Q, 
                      VT<double> &c)
  {
    const size_t d = A.rows();
    MT<double> dp(p.size(), p.size());
    genDiag(p, dp);

    MT<double> PN;
    PN.noalias() = dp * A.transpose();
    PN = A * PN;

    MT<double> M2;
    M2.noalias() = A * p;
    
    MT<double> M3;
    M3.noalias() = M2 * M2;

    MT<double> invert(PN.rows(), PN.cols());
    InvertLP(PN- M3, invert);
    Q = 1.0/d * invert;
    c.noalias() = A * p;

  }

  inline double KhachiyanAlgo(const MT<double> &A,
                       double eps,
                       size_t maxiter,
                       MT<double> &Q,
                       VT<double> &c)
  {
    VT<double> p(A.cols());
    p.setConstant(1.0/A.cols());

    MT<double> Ap;
    Lift(A, Ap);

    double ceps=eps*2;
    for (size_t i=0;  i<maxiter && ceps>eps; ++i)
    {
      ceps=KhachiyanIter(Ap, p);
    }

    KaInvertDual(A, p, Q, c);

    return ceps;


  }

  inline double KhachiyanAlgo(const std::set<MCPoint> &ss,
                       double eps,
                       size_t maxiter,
                       KhachiyanEllipsoid &res)
  {
    const size_t d=ss.begin()->p.size();
    MT<double> A(d, ss.size());

    size_t j=0;
    for (std::set<MCPoint>::const_iterator i=ss.begin();
         i != ss.end();
         ++i)
    {
      for(size_t k=0; k <d ;++k)
        A(k,j)=i->p[k];
      ++j;
    }

    MT<double> Q(d,d);
    VT<double> c(d);

    const double ceps=KhachiyanAlgo(A, eps, maxiter,
                                    Q, c);
    res.Q=Q;
    res.c=c;
    return ceps;
  }

#endif

//}
