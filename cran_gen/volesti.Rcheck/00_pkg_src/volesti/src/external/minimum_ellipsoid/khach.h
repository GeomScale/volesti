// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

//#include "khach1.h"
//#include "mcpoint1.h"
#include "mcpoint.h"
//#include "bnmin_main1.h"
//#include "bnmin_main2.h"

//#include "../bnmin_main.hxx"

//namespace Minim {

  namespace ublas=boost::numeric::ublas;

  struct KhachiyanEllipsoid
  {
      ublas::matrix<double> Q;
      ublas::vector<double> c;
  };

  template<class T>
  bool InvertMatrix(const ublas::matrix<T> &input,
                    ublas::matrix<T> &inverse)
  {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    matrix<T> A(input);
    pmatrix pm(A.size1());
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;
    inverse.assign(ublas::identity_matrix<T>(A.size1()));
    lu_substitute(A, pm, inverse);
    return true;
  }


  inline void InvertLP(const ublas::matrix<double> &Lambdap,
                ublas::matrix<double> &LpInv)
  {
    bool res=InvertMatrix(Lambdap, LpInv);
    if (not res)
    {
     // throw an error of your choice here!
     // throw MatrixErr("Could not invert matrix: ",
     //                 Lambdap);
    }
  }

  inline void Lift(const ublas::matrix<double> &A,
            ublas::matrix<double> &Ap)
  {
    Ap.resize(A.size1()+1,
              A.size2());
    ublas::matrix_range<ublas::matrix<double> >
      sub(Ap,
          ublas::range(0, A.size1()),
          ublas::range(0, A.size2()));
    sub.assign(A);
    ublas::row(Ap, Ap.size1()-1)=ublas::scalar_vector<double>(A.size2(),1.0);

  }

  inline void genDiag(const ublas::vector<double> &p,
               ublas::matrix<double> &res)
  {
    res.assign(ublas::zero_matrix<double>(p.size(),
                                          p.size()));
    for(size_t i=0; i<p.size(); ++i)
    {
      res(i,i)=p(i);
    }
  }

  inline void KaLambda(const ublas::matrix<double> &Ap,
                const ublas::vector<double> &p,
                ublas::matrix<double> &Lambdap)
  {

    ublas::matrix<double> dp(p.size(), p.size());
    genDiag(p, dp);

    dp=ublas::prod(dp, ublas::trans(Ap));
    Lambdap=ublas::prod(Ap,
                        dp);
  }

  inline double KhachiyanIter(const ublas::matrix<double> &Ap,
                       ublas::vector<double> &p)
  {
    /// Dimensionality of the problem
    const size_t d=Ap.size1()-1;

    ublas::matrix<double> Lp;
    ublas::matrix<double> M;
    KaLambda(Ap, p, Lp);
    ublas::matrix<double> ILp(Lp.size1(), Lp.size2());
    InvertLP(Lp, ILp);
    M=ublas::prod(ILp, Ap);
    M=ublas::prod(ublas::trans(Ap), M);

    double maxval=0;
    size_t maxi=0;
    for(size_t i=0; i<M.size1(); ++i)
    {
      if (M(i,i) > maxval)
      {
        maxval=M(i,i);
        maxi=i;
      }
    }
    const double step_size=(maxval -d - 1)/((d+1)*(maxval-1));
    ublas::vector<double> newp=p*(1-step_size);
    newp(maxi) += step_size;

    const double err= ublas::norm_2(newp-p);
    p=newp;
    return err;

  }

  inline void KaInvertDual(const ublas::matrix<double> &A,
                    const ublas::vector<double> &p,
                    ublas::matrix<double> &Q,
                    ublas::vector<double> &c
                    )
  {
    const size_t d=A.size1();
    ublas::matrix<double> dp(p.size(), p.size());
    genDiag(p, dp);

    ublas::matrix<double> PN=ublas::prod(dp, ublas::trans(A));
    PN=ublas::prod(A, PN);

    ublas::vector<double> M2=ublas::prod(A, p);
    ublas::matrix<double> M3=ublas::outer_prod(M2, M2);

    ublas::matrix<double> invert(PN.size1(), PN.size2());
    InvertLP(PN- M3, invert);

    Q.assign( 1.0/d *invert);
    c=ublas::prod(A, p);


  }

  inline double KhachiyanAlgo(const ublas::matrix<double> &A,
                       double eps,
                       size_t maxiter,
                       ublas::matrix<double> &Q,
                       ublas::vector<double> &c)
  {
    ublas::vector<double> p=ublas::scalar_vector<double>(A.size2(), 1.0)*(1.0/A.size2());

    ublas::matrix<double> Ap;
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
    ublas::matrix<double> A(d,
                            ss.size());

    size_t j=0;
    for (std::set<MCPoint>::const_iterator i=ss.begin();
         i != ss.end();
         ++i)
    {
      for(size_t k=0; k <d ;++k)
        A(k,j)=i->p[k];
      ++j;
    }

    ublas::matrix<double> Q(d,d);
    ublas::vector<double> c(d);

    const double ceps=KhachiyanAlgo(A, eps, maxiter,
                                    Q, c);
    res.Q=Q;
    res.c=c;
    return ceps;
  }

#endif

//}
