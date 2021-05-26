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
  using MTT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;

  template <class NT>
  using VTT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

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
  bool InvertMatrix(const MTT<T> &input,
                    MTT<T> &inverse)
  {
    inverse = input.inverse();
    return !is_nan(inverse);
  }


  inline void InvertLP(const MTT<double> &Lambdap,
                MTT<double> &LpInv)
  {
    bool res = InvertMatrix(Lambdap, LpInv);
    if (not res)
    {
     // throw an error of your choice here!
     // throw MatrixErr("Could not invert matrix: ",
     //                 Lambdap);
    }
  }

  inline void Lift(const MTT<double> &A, MTT<double> &Ap)
  {
    Ap.resize(A.rows()+1, A.cols());
    Ap.topLeftCorner(A.rows(), A.cols()) = A;
    Ap.row(Ap.rows()-1).setConstant(1.0); 
  }

  inline void genDiag(const VTT<double> &p, MTT<double> &res)
  {
    res.setZero(p.size(), p.size());

    for(size_t i=0; i<p.size(); ++i)
    {
      res(i,i)=p(i);
    }
  }

  inline void KaLambda(const MTT<double> &Ap,
                const VTT<double> &p,
                MTT<double> &Lambdap)
  {

    MTT<double> dp(p.size(), p.size());
    genDiag(p, dp);

    dp = dp * Ap.transpose();
    Lambdap.noalias() = Ap * dp;
  }

  inline double KhachiyanIter(const MTT<double> &Ap, VTT<double> &p)
  {
    /// Dimensionality of the problem
    const size_t d = Ap.rows()-1;

    MTT<double> Lp;
    MTT<double> M;
    KaLambda(Ap, p, Lp);
    MTT<double> ILp(Lp.rows(), Lp.cols());
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
    VTT<double> newp = p*(1-step_size);
    newp(maxi) += step_size;

    const double err= (newp-p).norm();
    p = newp;
    return err;

  }

  inline void KaInvertDual(const MTT<double> &A, 
                      const VTT<double> &p, 
                      MTT<double> &Q, 
                      VTT<double> &c)
  {
    const size_t d = A.rows();
    MTT<double> dp(p.size(), p.size());
    genDiag(p, dp);

    MTT<double> PN;
    PN.noalias() = dp * A.transpose();
    PN = A * PN;

    VTT<double> M2;
    M2.noalias() = A * p;
    
    MTT<double> M3;
    M3.noalias() = M2 * M2.transpose();

    MTT<double> invert(PN.rows(), PN.cols());
    InvertLP(PN- M3, invert);
    Q.noalias() = (invert/d);
    c.noalias() = A * p;

  }

  inline double KhachiyanAlgo(const MTT<double> &A,
                       double eps,
                       size_t maxiter,
                       MTT<double> &Q,
                       VTT<double> &c)
  {
    VTT<double> p(A.cols());
    p.setConstant(1.0/A.cols());

    MTT<double> Ap;
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
    MTT<double> A(d, ss.size());

    size_t j=0;
    for (std::set<MCPoint>::const_iterator i=ss.begin();
         i != ss.end();
         ++i)
    {
      for(size_t k=0; k <d ;++k)
        A(k,j)=i->p[k];
      ++j;
    }

    MTT<double> Q(d,d);
    VTT<double> c(d);

    const double ceps=KhachiyanAlgo(A, eps, maxiter,
                                    Q, c);
    res.Q=Q;
    res.c=c;
    return ceps;
  }

#endif

//}
