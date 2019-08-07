//
// Created by panagiotis on 5/8/2019.
//

#ifndef VOLESTI_TRUNCATED_EXPONENTIAL_H
#define VOLESTI_TRUNCATED_EXPONENTIAL_H

double expPDF(double x, double lambda) {
    return lambda*std::exp(-lambda*x);
}

double expCDF(double x, double lambda) {
    return 1 - std::exp(- lambda*x);
}

double expQuantile(double p, double lambda) {
    return -std::log(1 - p) / lambda;
}


/**
 * Produce a randomly generated number from a truncated distribution
 *
 * @param x
 * @param lambda
 * @param a lower bound
 * @param b upper bound
 * @return
 */
template <class RNGType>
double texp(double lambda, double a, double b, RNGType& rng) {
    boost::random::uniform_real_distribution<> urdist(0, 1);
    double u = urdist(rng);
    double cdfA = expCDF(a, lambda);
    double cdfB = expCDF(b, lambda);

    return expQuantile(cdfA + u*(cdfB - cdfA), lambda);
}

#endif //VOLESTI_TRUNCATED_EXPONENTIAL_H
