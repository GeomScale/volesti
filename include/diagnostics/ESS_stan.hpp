
#ifndef ESS_STAN_HPP
#define ESS_STAN_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <vector>
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

/*  The functions `get_good_size_2()`, `autocorrelation()` and  `autocovariance()`
    was taken from rstan code at https://github.com/stan-dev/rstan */

inline size_t get_good_size_2(size_t N) {
  // Find the optimal next size for the FFT so that
  // a minimum number of zeros are padded.
  
  if (N <= 2) {
    return(2);
  }
  size_t m;
  while (TRUE) {
    m = N;
    while ((m % 2) == 0){ 
      m = m / 2;
    }
    while ((m % 3) == 0){
       m = m / 3;
    }
    while ((m % 5) == 0) { 
      m = m / 5;
    }
    if (m <= 1) {
      return(N);
    }
    N = N + 1;
  }
}


/**
 * Write autocorrelation estimates for every lag for the specified
 * input sequence into the specified result using the specified FFT
 * engine. Normalizes lag-k autocorrelation estimators by N instead
 * of (N - k), yielding biased but more stable estimators as
 * discussed in Geyer (1992); see
 * https://projecteuclid.org/euclid.ss/1177011137. The return vector
 * will be resized to the same length as the input sequence with
 * lags given by array index.
 *
 * <p>The implementation involves a fast Fourier transform,
 * followed by a normalization, followed by an inverse transform.
 *
 * <p>An FFT engine can be created for reuse for type double with:
 *
 * <pre>
 *     Eigen::FFT<double> fft;
 * </pre>
 *
 * @tparam T Scalar type.
 * @param y Input sequence.
 * @param ac Autocorrelations.
 * @param fft FFT engine instance.
 */
template <typename T, typename DerivedA, typename DerivedB>
void autocorrelation(const Eigen::MatrixBase<DerivedA>& y,
                     Eigen::MatrixBase<DerivedB>& ac, Eigen::FFT<T>& fft) {
  size_t N = y.size();
  size_t M = get_good_size_2(N);
  size_t Mt2 = 2 * M;

  // centered_signal = y-mean(y) followed by N zeros
  Eigen::Matrix<T, Eigen::Dynamic, 1> centered_signal(Mt2);
  centered_signal.setZero();
  centered_signal.head(N) = y.array() - y.mean();

  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> freqvec(Mt2);
  fft.fwd(freqvec, centered_signal);
  // cwiseAbs2 == norm
  freqvec = freqvec.cwiseAbs2();

  Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> ac_tmp(Mt2);
  fft.inv(ac_tmp, freqvec);

  // use "biased" estimate as recommended by Geyer (1992)
  ac = ac_tmp.head(N).real().array() / (N * N * 2);
  ac /= ac(0);
}

/**
 * Write autocovariance estimates for every lag for the specified
 * input sequence into the specified result using the specified FFT
 * engine. Normalizes lag-k autocovariance estimators by N instead
 * of (N - k), yielding biased but more stable estimators as
 * discussed in Geyer (1992); see
 * https://projecteuclid.org/euclid.ss/1177011137. The return vector
 * will be resized to the same length as the input sequence with
 * lags given by array index.
 *
 * <p>The implementation involves a fast Fourier transform,
 * followed by a normalization, followed by an inverse transform.
 *
 * <p>This method is just a light wrapper around the three-argument
 * autocovariance function
 *
 * @tparam T Scalar type.
 * @param y Input sequence.
 * @param acov Autocovariances.
 */
template <typename T, typename DerivedA, typename DerivedB>
void autocovariance(const Eigen::MatrixBase<DerivedA>& y,
                    Eigen::MatrixBase<DerivedB>& acov) {
  Eigen::FFT<T> fft;
  autocorrelation(y, acov, fft);

  using boost::accumulators::accumulator_set;
  using boost::accumulators::stats;
  using boost::accumulators::tag::variance;

  accumulator_set<double, stats<variance>> acc;
  for (int n = 0; n < y.size(); ++n) {
    acc(y(n));
  }

  acov = acov.array() * boost::accumulators::variance(acc);
}



template <typename NT, typename VT, typename MT>
class ESSestimator {

private:
   unsigned int         num_draws, max_s, s, d, num_chains, jj; 
   VT                   cm_mean, cm_var, cv_mean, draws, var_plus, ess; 
   NT                   oldM, rho_hat_odd, rho_hat_even, mean_var, M2, delta, new_elem;
   MT                   acov_s_mean, rho_hat_s;
   Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, 1> acov;
   
public:
   ESSestimator() {}

   ESSestimator(unsigned int const& _ndraws, unsigned int const& _dim) 
   {
     num_draws = _ndraws;
     d = _dim;
     num_chains = 0;

     cm_mean.setZero(d);
     cm_var.setZero(d);
     cv_mean.setZero(d);
     var_plus.setZero(d);
     ess.setZero(d);
     draws.setZero(num_draws);
     acov_s_mean.setZero(num_draws-3, d);
     rho_hat_s.setZero(num_draws, d);
     acov = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, 1>(1);
   }

  void update_estimator(MT const& samples) 
  {
    num_chains++;
    for (int i = 0; i < d; i++)
    {
      draws = samples.row(i).transpose();
      autocovariance<double>(draws, acov(0));

      new_elem = draws.mean();
      delta = new_elem - cm_mean.coeff(i);
      cm_mean(i) += delta / NT(num_chains);
      cm_var(i) += delta * (new_elem - cm_mean(i));

      new_elem = acov(0)(0) * NT(num_draws) / (NT(num_draws) - 1.0);
      delta = new_elem - cv_mean.coeff(i);
      cv_mean(i) += delta / NT(num_chains);

      new_elem = acov(0)(1);
      delta = new_elem - acov_s_mean.coeff(0, i);
      acov_s_mean(0, i) += delta / NT(num_chains);
      jj = 1;
      while (jj < num_draws-4)
      {
        new_elem = acov(0)(jj+1);
        delta = new_elem - acov_s_mean.coeff(jj, i);
        acov_s_mean(jj, i) += delta / NT(num_chains);

        new_elem = acov(0)(jj+2);
        delta = new_elem - acov_s_mean.coeff(jj+1, i);
        acov_s_mean(jj+1, i) += delta / NT(num_chains);

        jj += 2;
      }
    }
    
  }

  void estimate_effective_sample_size()
  {
    rho_hat_s.setZero(num_draws, d); 

    var_plus = cv_mean * (NT(num_draws) - 1.0) / NT(num_draws);
    if (num_chains > 1) {
      VT cm_var_temp = cm_var * (1.0 / (NT(num_chains)-1.0));
      var_plus += cm_var_temp;
    }

    for (int i = 0; i < d; i++)
    {
      rho_hat_even = 1.0;
      rho_hat_s(0, i) = rho_hat_even;
      rho_hat_odd = 1 - (cv_mean.coeff(i) - acov_s_mean.coeff(0, i)) / var_plus.coeff(i);
      rho_hat_s(1, i) = rho_hat_odd;

      s = 1;
      while (s < (num_draws - 4) && (rho_hat_even + rho_hat_odd) > 0) {
        rho_hat_even = 1.0 - (cv_mean.coeff(i) - acov_s_mean.coeff(s, i)) / var_plus.coeff(i);
        rho_hat_odd = 1.0 - (cv_mean.coeff(i) - acov_s_mean.coeff(s+1, i)) / var_plus.coeff(i);
        if ((rho_hat_even + rho_hat_odd) >= 0) {
          rho_hat_s(s + 1, i) = rho_hat_even;
          rho_hat_s(s + 2, i) = rho_hat_odd;
        }
        s += 2;
      }

      max_s = s;
      // this is used in the improved estimate, which reduces variance
      // in antithetic case -- see tau_hat below
      if (rho_hat_even > 0) {
        rho_hat_s(max_s + 1, i) = rho_hat_even;
      }

      // Convert Geyer's initial positive sequence into an initial
      // monotone sequence
      for (jj = 1; jj <= max_s - 3; jj += 2) {
        if (rho_hat_s(jj + 1, i) + rho_hat_s.coeff(jj + 2, i) > rho_hat_s.coeff(jj - 1, i) + rho_hat_s.coeff(jj, i)) {
          rho_hat_s(jj + 1, i) = (rho_hat_s.coeff(jj - 1, i) + rho_hat_s.coeff(jj, i)) / 2.0;
          rho_hat_s(jj + 2, i) = rho_hat_s.coeff(jj + 1, i);
        }
      }

      NT num_total_draws = NT(num_chains) * NT(num_draws);
      NT tau_hat = -1.0 + 2.0 * rho_hat_s.col(i).head(max_s).sum() + rho_hat_s.coeff(max_s + 1, i);
      ess(i) = num_total_draws / tau_hat;
    }
    
  }

  VT get_effective_sample_size() 
  {
    return ess;
  }

};


#endif
