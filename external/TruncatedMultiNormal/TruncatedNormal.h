#pragma once
#ifndef MCQMC_INTEGRATION_TRUNCATED_NORMAL_H
#define MCQMC_INTEGRATION_TRUNCATED_NORMAL_H
/**
 * @file TruncatedNormal.h
 *
 * @brief Calculate Truncated Multivalual Normal Distribution
 *
 * @author Zdravko I. Botev (original R version)
 * @author Shinsuke Mori (port to C++)
 * @author Makoto Matsumoto
 * @author Mutsuo Saito (little bit modification for publish)
 *
 * Copyright (C) 2015 Zdravko I. Botev.
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 */

#include "DigitalNet2.h"
#include <vector>
#include <inttypes.h>

namespace MCQMCIntegration {

    /**
     * structure of the result of truncatedNormal.
     */
    struct TruncatedNormalResult {
        /** estimated value of probability Pr(@c lower < X < @c upper). */
        double probability;
        /** absolute error. */
        double absoluteError;
        /** relative error. */
        double relativeError;
        /** theoretical upper bound on true Pr(@c lower < X < @c upper). */
        double upperBound;
        /** calculation success or not. */
        bool success;
    };

    /**
     * compute Multivariate Normal Distribution by Quasi Monte-Carlo method.
     *
     * computes an estimator and a deterministic upper bound of the
     * probability P r(@c lower < X < @c upper), where X is a
     * zero-mean multivariate normal vector with covariance matrix Σ,
     * that is, X is drawn from N(0, Σ) infinite values for vectors @c
     * upper and @c lower are accepted.
     *
     * @note parameter @c dn, digital net should have dimension s =
     * lower.size - 1. And s should be an integer 4 <= s <= 10,
     * therefore, lower.size size should be 5 <= size <= 11.
     * Moreover, The F2-dimension of element of digitalnet, @c m
     * should be 10 <= @c m <= 18.
     *
     * @note @c probability should be one of {0.95, 0.99, 0.999, 0.9999}
     *
     * @param[in] lower lower truncation limit.
     * @param[in] upper upper truncation limit.
     * @param[in] sigma covariance matrix of N(0, Σ)
     * @param[in] number Monte Carlo simulation effort — the larger the n,
     * the smaller the relative error of the estimator.
     * @param[in,out] dn digital net.
     * @param[in] probability expected probability such that the
     * result is in the range [p-σ, p+σ].
     * @return estimated value of probability and other information.
     */
    TruncatedNormalResult truncatedNormal(
        const std::vector<double>& lower,
        const std::vector<double>& upper,
        const std::vector< std::vector<double> >& sigma,
        uint64_t number,
        DigitalNet<uint64_t>& dn,
        double probability);

    /**
     * compute Multivariate Normal Distribution by Quasi Monte-Carlo method.
     *
     * computes an estimator and a deterministic upper bound of the
     * probability P r(@c lower < X < @c upper), where X is a
     * zero-mean multivariate normal vector with covariance matrix Σ,
     * that is, X is drawn from N(0, Σ) infinite values for vectors @c
     * upper and @c lower are accepted.
     *
     * @c dnid is one of
     * @li @c DigitalNetID::NXLW Niederreiter-Xing Low WAFOM
     * @li @c DigitalNetID::SOLW Sobol Low WAFOM
     *
     * @note The F2-dimension of element of digitalnet, @c m
     * should be 10 <= @c m <= 18.
     *
     * @c probability should be one of {0.95, 0.99, 0.999, 0.9999}
     *
     * @param[in] lower lower truncation limit.
     * @param[in] upper upper truncation limit.
     * @param[in] sigma covariance matrix of N(0, Σ)
     * @param[in] number Monte Carlo simulation effort — the larger
     * the @c number, the smaller the relative error of the estimator.
     * @param[in] dnid digital netid.
     * @param[in] m F2 dimension of element of digital net.
     * @param[in] probability expected probability such that the
     * result is in the range [p - @c absoluteError, p + @c absoluteError].
     * @return estimated value of probability and other information.
     */
    TruncatedNormalResult truncatedNormal(
        const std::vector<double>& lower,
        const std::vector<double>& upper,
        const std::vector< std::vector<double> >& sigma,
        uint64_t number,
        DigitalNetID dnid,
        uint32_t m,
        double probability);

    /**
     * compute Multivariate Normal Distribution by Monte-Carlo method.
     *
     * computes an estimator and a deterministic upper bound of the
     * probability P r(@c lower < X < @c upper), where X is a
     * zero-mean multivariate normal vector with covariance matrix Σ,
     * that is, X is drawn from N(0, Σ) infinite values for vectors @c
     * upper and @c lower are accepted; Monte Carlo method uses sample
     * size @c number;
     *
     * @c probability should be one of {0.95, 0.99, 0.999, 0.9999}
     *
     * @param[in] lower lower truncation limit.
     * @param[in] upper upper truncation limit.
     * @param[in] sigma covariance matrix of N(0, Σ)
     * @param[in] trialNumber Monte Carlo simulation effort — the larger
     * the @c number, the smaller the relative error of the estimator.
     * @param[in] sampleNumber sample number per a trial.
     * @param[in] probability expected probability such that the
     * result is in the range [p - @c absoluteError, p + @c absoluteError].
     * @return estimated value of probability and other information.
     */
    TruncatedNormalResult truncatedNormalMC(
        const std::vector<double>& lower,
        const std::vector<double>& upper,
        const std::vector< std::vector<double> >& sigma,
        uint32_t trialNumber,
        uint32_t sampleNumber,
        double probability);

}
#endif // MCQMC_INTEGRATION_TRUNCATED_NORMAL_H
