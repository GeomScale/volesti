// Code from https://github.com/rzrsk/vaidya-walk

// Modified by Apostolos Chalkis, as part of Google Summer of Code 2019 program.


#ifndef PWALK_UTIL_MATH_FUNCTIONS_HPP_
#define PWALK_UTIL_MATH_FUNCTIONS_HPP_

#include <Eigen/Dense>
#include <boost/random.hpp>
#include <cmath>

//namespace pwalk {

// Define random number generator type
typedef boost::mt19937 rng_t;

template <typename Dtype>
void sample_gaussian(const int n, const Dtype a,
                     const Dtype sigma, Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& r) {
    static rng_t gen(1234567);
    static boost::normal_distribution<Dtype> random_distribution(a, sigma);
    static boost::variate_generator<rng_t&, boost::normal_distribution<Dtype> >
    variate_generator(gen, random_distribution);

    for (int i = 0; i < n; ++i) {
        r[i] = variate_generator();
    }
}

template <typename Dtype>
Dtype rng_uniform(const Dtype a, const Dtype b) {
    static rng_t gen(1234567);
    static boost::uniform_real<Dtype> random_distribution(a, b);
    static boost::variate_generator<rng_t&, boost::uniform_real<Dtype> >
    variate_generator(gen, random_distribution);

    return variate_generator();
}

template <typename Dtype>
Dtype gaussian_density(const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& x, const Eigen::Matrix<Dtype, Eigen::Dynamic, 1>& mu, const Eigen::Matrix<Dtype, Eigen::Dynamic, Eigen::Dynamic>& sqrt_cov) {
    Eigen::Matrix<Dtype, Eigen::Dynamic, 1> c = sqrt_cov * (x - mu);
    return std::exp(-0.5*c.dot(c)) * sqrt_cov.determinant();
}


#endif // PWALK_UTIL_MATH_FUNTIONS_HPP_
