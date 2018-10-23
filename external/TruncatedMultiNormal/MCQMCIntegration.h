#pragma once
#ifndef MCQMC_INTEGRATION_HPP
#define MCQMC_INTEGRATION_HPP
/**
 * @file MCQMCIntegration.h
 *
 * @brief Monte-Carlo and Quasi Monte-Carlo Integration
 *
 * @author Shinsuke Mori (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Mutsuo Saito
 *
 * Copyright (C) 2017 Shinsuke Mori, Makoto Matsumoto, Mutsuo Saito
 * and Hiroshima University.
 * All rights reserved.
 *
 * The GPL ver.3 is applied to this software, see
 * COPYING
 */

#include "DigitalNet2.h"
#include <random>

namespace MCQMCIntegration {

    /*
     * calculate variance
     *
     * See "Online algorithm" at
     * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
     */
    class OnlineVariance {
    public:
        OnlineVariance();
        void addData(const double x) {
            n++;
            df = n - 1;
            double delta = x - mean;
            mean += delta / static_cast<double>(n);
            M2 += delta * (x - mean);
        }
        double getMean() const;
        double unbiasedVar() const;
        double var() const;
        double absErr(const int prob) const;
        double relErr(const int prob) const;
    private:
        int n;
        int df;
        double mean;
        double M2;
    };

    /**
     * Result Structure of Numeric Integration
     *
     */
    struct MCQMCResult {
        double value;
        double error;
    };

    /*
     * Monte-Carlo Integration.
     *
     * @tperm I integrand function class.
     * @tparm R Random number generator class.
     * @tparm D Distribution class for Random number generator.
     *
     * @param[in] s dimension of integration area @b R.
     * @param[in] m sample number per a trial.
     * @param[in] N number of trials.
     * @param[in,out] integrand integrand function class, which should have
     * double operator()(double[]).
     * @param[in,out] rand random number generator.
     * @param[in,out] dist random number distribution class.
     * @param[in] probability expected probability of returned value x is
     * between x - absolute error and x + absolute error. this should be
     * one of {95, 99, 999, 9999}.
     * @return MCQMCResult.
     */
    template<typename I, typename R, typename D>
        MCQMCResult monte_carlo_integration(uint32_t s,
                                            uint32_t m,
                                            uint32_t N,
                                            I& integrand,
                                            R& rand,
                                            D& dist,
                                            int probability = 99)
    {
        double point[s];
        OnlineVariance eachintval;
        uint32_t cnt = 0;
        do {
            OnlineVariance intsum;
            for (uint32_t j = 0; j < m; ++j) {
                for (uint32_t i = 0; i < s; ++i) {
                    point[i] = dist(rand);
                }
                intsum.addData(integrand(point));
            }
            eachintval.addData(intsum.getMean());
            ++cnt;
        } while ( cnt < N );
        return MCQMCResult({eachintval.getMean(),
                    eachintval.absErr(probability)});
    }

    /*
     * Quasi Monte-Carlo Integration
     *
     * @tperm I integrand function class
     * @tparm D DigitalNet class for Quasi Monete-Carlo integration.
     *
     * @param[in] N number of trials.
     * @param[in,out] integrand integrand function class, which should have
     * double operator()(double[]).
     * @param[in,out] digitalNet digital net class.
     * @param[in] probability expected probability of returned value x is
     * between x - absolute error and x + absolute error. this should be
     * one of {95, 99, 999, 9999}.
     * @return MCQMCResult.
     */
    template<typename I, typename D>
        MCQMCResult quasi_monte_carlo_integration(uint32_t N,
                                                  I& integrand,
                                                  D& digitalNet,
                                                  int probability = 99)
    {
        uint32_t m = digitalNet.getM();
        digitalNet.pointInitialize();
        OnlineVariance eachintval;
        uint32_t cnt = 0;
        do {
            OnlineVariance intsum;
            uint64_t max = 1;
            max = max << m;
            for (uint64_t j = 0; j < max; ++j) {
                intsum.addData(integrand(digitalNet.getPoint()));
                digitalNet.nextPoint();
            }
            eachintval.addData(intsum.getMean());
            cnt++;
        } while ( cnt < N );
        return MCQMCResult({eachintval.getMean(),
                    eachintval.absErr(probability)});
    }

    /*
     * Quasi Monte-Carlo Integration
     *
     * @tperm I integrand function class
     *
     * @param[in] N number of trials.
     * @param[in,out] integrand integrand function class, which should have
     * double operator()(double[]).
     * @param[in,out] digitalNet digital net class.
     * @param[in] probability expected probability of returned value x is
     * between x - absolute error and x + absolute error. this should be
     * one of {95, 99, 999, 9999}.
     * @return MCQMCResult.
     */
    template<typename I>
        MCQMCResult quasi_monte_carlo_integration(uint32_t N,
                                                  I& integrand,
                                                  DigitalNetID digitalNetId,
                                                  uint32_t s,
                                                  uint32_t m,
                                                  int probability)
    {
        DigitalNet<uint64_t> digitalNet(digitalNetId, s, m);
        digitalNet.pointInitialize();
        OnlineVariance eachintval;
        uint32_t cnt = 0;
        do {
            OnlineVariance intsum;
            uint64_t max = 1;
            max = max << m;
            for (uint64_t j = 0; j < max; ++j) {
                intsum.addData(integrand(digitalNet.getPoint()));
                digitalNet.nextPoint();
            }
            eachintval.addData(intsum.getMean());
            cnt++;
        } while ( cnt < N );
        return MCQMCResult({eachintval.getMean(),
                    eachintval.absErr(probability)});
    }
}
#endif // MCQMC_INTEGRATION_HPP
