// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef ESS_UPDATER_UTILS_HPP
#define ESS_UPDATER_UTILS_HPP

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
  while (true) {
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



#endif