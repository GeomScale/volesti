/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2009

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file mcpoint.hxx

   Structure representing a Markov Chain point
*/

#ifndef _BNMIN1_MCPOINT_HXX__
#define _BNMIN1_MCPOINT_HXX__

#include <vector>
#include <list>
#include <set>

namespace Minim {

  /** \brief Data to be recorded at each point in an MCMC distribution
  */
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
    MCPoint(void);

    /** \Construct with supplied position vector
     */
    MCPoint(const std::vector<double> &p);

    /** \brief The parameter vector has n values
     */
    MCPoint(size_t np);

    MCPoint(const MCPoint &other);

    MCPoint & operator=(const MCPoint &other);
    
    
  };

  /** \brief Define ordering of MCPoints on basis of their likelihood

      Should be useful for nested sampling
  */
  inline bool operator< (const MCPoint &a, const MCPoint &b)
  {
    return a.ll < b.ll;
  }

  /** \brief Weighted posterior point

      As produced by non markov chain posterior exploration
   */
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

  /** \brief Calculate the first moment of each parameter
   */
  void moment1(const std::list<WPPoint> &l,
	       std::vector<double> &res);

  void moment1(const std::list<WPPoint> &l,
	       double Z,
	       std::vector<double> &res);

  /** \brief Calculate the second moment of each parameter
   */
  void moment2(const std::list<WPPoint> &l,
	       const std::vector<double> &m1,
	       std::vector<double> &res);

  void moment2(const std::list<WPPoint> &l,
	       const std::vector<double> &m1,
	       double Z,
	       std::vector<double> &res);


  void moment1(const std::set<MCPoint> &s,
	       std::vector<double> &res);

  void moment2(const std::set<MCPoint> &s,
	       const std::vector<double> &m1,
	       std::vector<double> &res);

  /** \brief Outer moment 2, i.e., covariances
   */
  void omoment2(const std::set<MCPoint> &s,
		const std::vector<double> &m1,
		std::vector<double> &res);

  /**  Covariances, leaving the m1 to be calculated internally
   */
  void omoment2(const std::set<MCPoint> &s,
		std::vector<double> &res);

  /** Standard eviation of each coordiante separately
   */
  void StdDev(const std::set<MCPoint> &s,
	      std::vector<double> &res);
  

  /** \brief Compute the principal directions from the covariance
      matrix

      \param eigvals The eigen values
      \param eigvects The eigen vectors
      
   */
  void principalCV(const std::vector<double> &cv,
		   std::vector<double> &eigvals,
		   std::vector<double> &eigvects);


  /** \brief Compute the histogram of the posterior from weighted
      sample points

      \param low The low boundary of the region to be binned in each
      dimension
      
      \param high The high boundary of the region to be binned 

      \param nbins Number of bins in each direction

      \param res Results are stored here
      
   */
  void postHist(const std::list<WPPoint> &l,
		double Z,
		const std::vector<double> &low,
		const std::vector<double> &high,
		size_t nbins,
		std::vector<double> &res);

  /** Calculate the marginalised probability histogram

      \param pi Parameter index to marginaise to

      \param Z Evidence value, used to normalize the histogram. If
      this is set to one, the integral will be the evidence rather
      than the marginalised probability given the hypothesis

      \param high The high boundary of region to histogram, i.e., the
      *top* of the highest bin
      
   */
  void marginHist(const std::list<WPPoint> &l,
		  size_t pi,
		  double Z,
		  double low,
		  double high,
		  size_t nbins,
		  std::vector<double> &res);

  /** Calculate the marginalisation of the posterior to two dimensions

      \param Z Evidence value, used to normalize the histogram. If
      this is set to one, the integral will be the evidence rather
      than the marginalised probability given the hypothesis

   */
  void marginHist2D(const std::list<WPPoint> &l,
		    double Z,
		    size_t i,
		    double ilow,
		    double ihigh,
		    size_t j,
		    double jlow,
		    double jhigh,
		    size_t nbins,
		    std::vector<double> &res);


}

#endif
