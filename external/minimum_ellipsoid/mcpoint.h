// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// This file is converted from BNMin1 (https://www.mrao.cam.ac.uk/~bn204/oof/bnmin1.html) by Apostolos Chalkis

// Original copyright notice:

/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2009

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file mcpoint.cxx
*/
#ifndef MCPOINT_H
#define MCPOINT_H

#include <vector>
#include <list>
#include <set>

#include <cmath>
#include <algorithm>

//exclude gsl library  Apostolos Chalkis
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_eigen.h>

//#include "mcpoint1.h"
//#include "mcpoint2.h"
#include "bnmin_main.h"
//#include "bnmin_main2.h"

//namespace Minim {
    struct MCPoint
    {
        /// The actual parameters
        std::vector<double> p;
        /// Log-likelihood of this point
        double ll;
        /// A vector to store derived quantities at sample of the
        /// distribution
        std::vector<double> fval;

        /// Default constructor allowed, fill in the data later
        MCPoint(void):
                p(0),
                ll(-9999),
                fval(0)
        {
        }

        /** \Construct with supplied position vector
         */
        MCPoint(const std::vector<double> &p):
                p(p),
                ll(-9999),
                fval(0)
        {
        }

        /** \brief The parameter vector has n values
         */
        MCPoint(size_t np):
                p(np),
                ll(-9999),
                fval(0)
        {
        }

        MCPoint(const MCPoint &other):
                p(other.p),
                ll(other.ll),
                fval(other.fval)
        {
        }

        MCPoint & operator=(const MCPoint &other)
        {
            p=other.p;
            ll=other.ll;
            fval=other.fval;
            return *this;
        }


    };

    inline bool operator< (const MCPoint &a, const MCPoint &b)
    {
        return a.ll < b.ll;
    }

    struct WPPoint:
            public MCPoint
    {
        /** \brief Weight factor
         */
        double w;

        WPPoint(void):
                w(0.0)
        {
        }

        WPPoint(const std::vector<double> &p,
                double w):
                MCPoint(p),
                w(w)
        {
        }

        /** \brief Construct from MCPoint and a supplied weight
         */
        WPPoint(const MCPoint &mp,
                double w):
                MCPoint(mp),
                w(w)
        {
        }

    };

    /*
  MCPoint::MCPoint(void):
    p(0),
    ll(-9999),
    fval(0)
  {
  }
  
  MCPoint::MCPoint(const std::vector<double> &p):
    p(p),
    ll(-9999),
    fval(0)
  {
  }  

  MCPoint::MCPoint(size_t np):
    p(np),
    ll(-9999),
    fval(0)
  {
  }

  MCPoint::MCPoint(const MCPoint &other):
    p(other.p),
    ll(other.ll),
    fval(other.fval)
  {
  }

  MCPoint &MCPoint::operator=(const MCPoint &other)
  {
    p=other.p;
    ll=other.ll;
    fval=other.fval;
    return *this;
  }*/


  inline void moment1(const std::list<WPPoint> &l,
	       std::vector<double > &res)
  {
    const size_t n=l.begin()->p.size();
    res=std::vector<double>(n, 0.0);
    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	  res[j]+= (i->p[j] * i->w * exp(- i->ll));
      }
    }
  }

  inline void moment1(const std::list<WPPoint> &l,
	       double Z,
	       std::vector<double> &res)
  {
    moment1(l,res);
    for(size_t j=0; j<res.size(); ++j)
      res[j] /= Z;
  }


  inline void moment2(const std::list<WPPoint> &l,
	       const std::vector<double > &m1,
	       std::vector<double > &res)
  {
    const size_t n=m1.size();
    res=std::vector<double>(n, 0.0);
    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	res[j]+= ( pow(i->p[j]-m1[j],2.0) * i->w * exp(- i->ll));
      }
    }
  }

  inline void moment2(const std::list<WPPoint> &l,
	       const std::vector<double> &m1,
	       double Z,
	       std::vector<double> &res)
  {
    moment2(l, m1, res);
    for(size_t j=0; j<res.size(); ++j)
      res[j] /= Z;    
  }

  inline void moment1(const std::set<MCPoint> &s,
	       std::vector<double> &res)
  {
    const size_t n=s.begin()->p.size();
    res=std::vector<double>(n, 0.0);

    size_t N=0;
    for(std::set<MCPoint>::const_iterator i=s.begin();
	i!= s.end();
	++i)
    {
      if(i->p.size() != n)
      {
	throw NParsErr("moment1", n, i->p.size());
      }
      for (size_t j=0; j<n; ++j)
      {
	res[j]+= (i->p[j]);
      }
      ++N;
    }
    
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]/=N;
    }
  }

  inline void moment2(const std::set<MCPoint> &s,
	       const std::vector<double> &m1,
	       std::vector<double> &res)
  {
    const size_t n=m1.size();
    res=std::vector<double>(n, 0.0);

    size_t N=0;
    for(std::set<MCPoint>::const_iterator i=s.begin();
	i!= s.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	res[j]+= pow(i->p[j]-m1[j], 2);
      }
      ++N;
    }
    
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]/=N;
    }
  }


  inline void omoment2(const std::set<MCPoint> &s,
		const std::vector<double > &m1,
		std::vector<double > &res)
  {
    const size_t n=m1.size();
    res=std::vector<double>(n*n, 0.0);

    size_t N=0;
    for(std::set<MCPoint>::const_iterator i=s.begin();
	i!= s.end();
	++i)
    {
      for (size_t j=0; j<n; ++j)
      {
	for(size_t k=0; k<n; ++k)
	{
	  res[j*n+k] += (i->p[j]-m1[j])*(i->p[k]-m1[k]);
	}
      }
      ++N;
    }
    
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]/=N;
    }

  }

  inline void omoment2(const std::set<MCPoint> &s,
		std::vector<double> &res)
  {
    std::vector<double> m1;
    moment1(s, m1);
    omoment2(s, m1, res);
  }


  inline void StdDev(const std::set<MCPoint> &s,
	      std::vector<double > &res)
  {
    std::vector<double> m1, m2;
    moment1(s, m1);
    moment2(s, m1, m2);
    res.resize(m2.size());
    for(size_t j=0; j<res.size(); ++j)
    {
      res[j]=pow(m2[j],0.5);
    }
  }

  //exlude principalCV function. Not needed + exclude gsl lib. Apostolos Chalkis
/*  void principalCV(const std::vector<double> &cv,
		   std::vector<double> &eigvals,
		   std::vector<double> &eigvects)
  {
    const size_t n=sqrt(cv.size());
    gsl_matrix_view m
      = gsl_matrix_view_array (const_cast<double*>(&cv[0]), n, n);

    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);

    gsl_eigen_symmv_workspace * w =
      gsl_eigen_symmv_alloc (n);

    gsl_eigen_symmv (&m.matrix,
		     eval,
		     evec,
		     w);

    gsl_eigen_symmv_free (w);
    
    gsl_eigen_symmv_sort (eval, 
			  evec,
			  GSL_EIGEN_SORT_ABS_ASC);

    eigvals.resize(n);
    eigvects.resize(n*n);
    for(size_t j=0; j<n; ++j)
    {
      eigvals[j]=gsl_vector_get (eval, j);
      for(size_t i=0; i<n; ++i)
      {
	eigvects[j*n+i]= gsl_matrix_get(evec, i,j);
      }
    }
    
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
  }*/

  inline void postHist(const std::list<WPPoint> &l,
		double Z,
		const std::vector<double> &low,
		const std::vector<double> &high,
		size_t nbins,
		std::vector<double> &res)
  {
    const size_t ndim=low.size();

    //res.resize(pow(nbins, static_cast<size_t>(ndim)));
    res.resize( static_cast<int>( pow(static_cast<long double>(nbins), static_cast<long double>(ndim)) ) );
    std::fill(res.begin(), res.end(), 0.0);
    

    std::vector<double> deltas(ndim);
    for(size_t i=0; i<ndim; ++i)
    {
      deltas[i]=(high[i]-low[i])/nbins;
    }

    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      bool inside=true;
      size_t k=0;
      for (size_t j=0; j<ndim; ++j)
      {
	int dimi = int((i->p[j]-low[j])/deltas[j]);
	if (dimi >= 0 and dimi < (int)nbins)
	{
	  k+= dimi * static_cast<int>( pow(static_cast<long double>(nbins), static_cast<long double>(ndim-j-1)) );
	}
	else
	{
	  inside=false;
	}
      }
      if (inside)
      {
	res[k]+= i->w * exp(- i->ll);
      }
    }
  }


  inline void marginHist(const std::list<WPPoint> &l,
		  size_t pi,
          double Z,
          double low,
          double high,
		  size_t nbins,
		  std::vector<double > &res)
  {
    res.resize(nbins);
    std::fill(res.begin(), res.end(), 
	      0.0);

    const double d=(high-low)/nbins;
    for(std::list<WPPoint>::const_iterator i=l.begin();
	i!= l.end();
	++i)
    {
      int k=int((i->p[pi]-low)/d);
      if (k > 0 and k < (int)nbins)
      {
	res[k]+= i->w * exp(- i->ll);
      }
    }

    for(size_t i=0; i<res.size(); ++i)
    {
      res[i]/=Z;
    }
  }

  inline void marginHist2D(const std::list<WPPoint> &l,
		    double Z,
		    size_t i,
		    double ilow,
		    double ihigh,
		    size_t j,
		    double jlow,
		    double jhigh,
		    size_t nbins,
		    std::vector<double> &res)
  {
    // Two dimensions only
    res.resize( static_cast<int>( pow(static_cast<long double>(nbins), static_cast<long double>(2)) ) );
    std::fill(res.begin(), res.end(), 
	      0.0);
    const double idelta=(ihigh-ilow)/nbins;
    const double jdelta=(jhigh-jlow)/nbins;
    
    for(std::list<WPPoint>::const_iterator p=l.begin();
	p!= l.end();
	++p)
    {
      
      int dimi = int((p->p[i]-ilow)/idelta);
      int dimj = int((p->p[j]-jlow)/jdelta);
      
      if (dimi >= 0 and   dimi<((int)nbins)  and   dimj >= 0 and  dimj < ((int)nbins))
      {
	const size_t k= dimi*nbins + dimj;
	res[k]+= p->w * exp(- p->ll);
      }
      
    }
  }
//}

#endif
