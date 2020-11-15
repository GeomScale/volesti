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

size_t fft_next_good_size(size_t N) {
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
  size_t M = fft_next_good_size(N);
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

/**
 * Computes the effective sample size (ESS) for the specified
 * parameter across all kept samples.  The value returned is the
 * minimum of ESS and the number_total_draws *
 * log10(number_total_draws).
 *
 * See more details in Stan reference manual section "Effective
 * Sample Size". http://mc-stan.org/users/documentation
 *
 * Current implementation assumes draws are stored in contiguous
 * blocks of memory.  Chains are trimmed from the back to match the
 * length of the shortest chain.  Note that the effective sample size
 * can not be estimated with less than four draws.
 *
 * @param draws stores pointers to arrays of chains
 * @param sizes stores sizes of chains
 * @return effective sample size for the specified parameter
 */
 template <typename VT>
inline double compute_effective_sample_size(VT draws){
                                            // sizes) {
  int num_chains = 1;
  size_t num_draws = draws.size();
  //for (int chain = 1; chain < num_chains; ++chain) {
  //  num_draws = std::min(num_draws, sizes[chain]);
  //}

  //if (num_draws < 4) {
  //  return std::numeric_limits<double>::quiet_NaN();
  //}

  // check if chains are constant; all equal to first draw's value
  //bool are_all_const = false;
  //Eigen::VectorXd init_draw = Eigen::VectorXd::Zero(num_chains);
  /*
  for (int chain_idx = 0; chain_idx < num_chains; chain_idx++) {
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> draw(
        draws[chain_idx], sizes[chain_idx]);

    for (int n = 0; n < num_draws; n++) {
      if (!boost::math::isfinite(draw(n))) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    init_draw(chain_idx) = draw(0);

    if (draw.isApproxToConstant(draw(0))) {
      are_all_const |= true;
    }
  }

  if (are_all_const) {
    // If all chains are constant then return NaN
    // if they all equal the same constant value
    if (init_draw.isApproxToConstant(init_draw(0))) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }*/

  Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, 1> acov(num_chains);
  Eigen::VectorXd chain_mean(num_chains);
  Eigen::VectorXd chain_var(num_chains);
  //for (int chain = 0; chain < num_chains; ++chain) {
  //Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>> draw(
        //draws[chain], sizes[chain]);
  int chain = 0;
  autocovariance<double>(draws, acov(chain));
  chain_mean(chain) = draws.mean();
  chain_var(chain) = acov(chain)(0) * num_draws / (num_draws - 1);
  //}

  double mean_var = chain_var.mean();
  double var_plus = mean_var * (num_draws - 1) / num_draws;
  //if (num_chains > 1)
   // var_plus += math::variance(chain_mean);
  Eigen::VectorXd rho_hat_s(num_draws);
  rho_hat_s.setZero();
  Eigen::VectorXd acov_s(num_chains);
  for (int chain = 0; chain < num_chains; ++chain)
    acov_s(chain) = acov(chain)(1);
  double rho_hat_even = 1.0;
  rho_hat_s(0) = rho_hat_even;
  double rho_hat_odd = 1 - (mean_var - acov_s.mean()) / var_plus;
  rho_hat_s(1) = rho_hat_odd;

  // Convert raw autocovariance estimators into Geyer's initial
  // positive sequence. Loop only until num_draws - 4 to
  // leave the last pair of autocorrelations as a bias term that
  // reduces variance in the case of antithetical chains.
  size_t s = 1;
  while (s < (num_draws - 4) && (rho_hat_even + rho_hat_odd) > 0) {
    for (int chain = 0; chain < num_chains; ++chain)
      acov_s(chain) = acov(chain)(s + 1);
    rho_hat_even = 1 - (mean_var - acov_s.mean()) / var_plus;
    for (int chain = 0; chain < num_chains; ++chain)
      acov_s(chain) = acov(chain)(s + 2);
    rho_hat_odd = 1 - (mean_var - acov_s.mean()) / var_plus;
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_s(s + 1) = rho_hat_even;
      rho_hat_s(s + 2) = rho_hat_odd;
    }
    s += 2;
  }

  int max_s = s;
  // this is used in the improved estimate, which reduces variance
  // in antithetic case -- see tau_hat below
  if (rho_hat_even > 0)
    rho_hat_s(max_s + 1) = rho_hat_even;

  // Convert Geyer's initial positive sequence into an initial
  // monotone sequence
  for (int s = 1; s <= max_s - 3; s += 2) {
    if (rho_hat_s(s + 1) + rho_hat_s(s + 2) > rho_hat_s(s - 1) + rho_hat_s(s)) {
      rho_hat_s(s + 1) = (rho_hat_s(s - 1) + rho_hat_s(s)) / 2;
      rho_hat_s(s + 2) = rho_hat_s(s + 1);
    }
  }

  double num_total_draws = num_chains * num_draws;
  // Geyer's truncated estimator for the asymptotic variance
  // Improved estimate reduces variance in antithetic case
  double tau_hat = -1 + 2 * rho_hat_s.head(max_s).sum() + rho_hat_s(max_s + 1);
  return std::min(num_total_draws / tau_hat,
                  num_total_draws * std::log10(num_total_draws));
}





#endif
