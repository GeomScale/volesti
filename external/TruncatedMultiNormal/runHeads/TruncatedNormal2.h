/**
 * @file TruncatedNormal.cpp
 *
 * @brief Calculate Truncated Multivalual Normal Distribution
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

#ifndef TRUNCATEDNORMALTWO_H
#define TRUNCATEDNORMALTWO_H

//#define DEBUG 1

#include <TruncatedNormal.h>
#include <MCQMCIntegration2.h>
#include <BotevAlgorithm2.h>
#include <random>

namespace {
    int probToInt(double probability) {
        double x = 1.0 - probability;
        if (x < 0.00011) {
            return 9999;
        } else if (x < 0.0011) {
            return 999;
        } else if (x < 0.011) {
            return 99;
        } else {
            return 95;
        }
    }
}

namespace MCQMCIntegration {
    using namespace std;

    TruncatedNormalResult truncatedNormal(
        const std::vector<double>& lower,
        const std::vector<double>& upper,
        const std::vector< std::vector<double> >& sigma,
        uint64_t number,
        DigitalNet<uint64_t>& dn,
        double probability)
    {
#if defined(DEBUG)
        cout << "truncatedNormal start" << endl;
#endif
        TruncatedNormalResult result;
        uint32_t size = lower.size();
        if (size != upper.size() || size != sigma.size() ||
            size != sigma[0].size()) {
            cerr << "size mismatch" << endl;
            result.success = false;
            return result;
        }
        if (size != dn.getS() + 1) {
            cerr << "digital net size mismatch" << endl;
            result.success = false;
            return result;
        }
#if defined(DEBUG)
        cout << "before ba constructor" << endl;
#endif
        BotevAlgorithm ba(lower, upper, sigma);
#if defined(DEBUG)
        cout << "after ba constructor" << endl;
#endif
        if (!ba.checkAndPrepare()) {
            result.success = false;
            return result;
        }
#if defined(DEBUG)
        cout << "after ba checkAndPrepare" << endl;
#endif
        MCQMCResult qmcresult;
        int p = probToInt(probability);
        qmcresult = quasi_monte_carlo_integration(number, ba, dn, p);
        result.probability = qmcresult.value;
        result.absoluteError = qmcresult.error;
        result.relativeError = result.absoluteError / result.probability;
        result.success = true;
        result.upperBound = ba.getUpperBound();
        return result;
    }

    TruncatedNormalResult truncatedNormal(
        const std::vector<double>& lower,
        const std::vector<double>& upper,
        const std::vector< std::vector<double> >& sigma,
        uint64_t number,
        DigitalNetID dnid,
        uint32_t m,
        double probability)
    {
        uint32_t s = lower.size();
        DigitalNet<uint64_t> dn(dnid, s - 1, m);
        return truncatedNormal(lower, upper, sigma, number, dn, probability);
    }

    TruncatedNormalResult truncatedNormalMC(
        const std::vector<double>& lower,
        const std::vector<double>& upper,
        const std::vector< std::vector<double> >& sigma,
        uint32_t trialNumber,
        uint32_t sampleNumber,
        double probability)
    {
        TruncatedNormalResult result;
        BotevAlgorithm ba(lower, upper, sigma);
        if (!ba.checkAndPrepare()) {
            result.success = false;
            return result;
        }
        MCQMCResult mcresult;
        int p = probToInt(probability);
        uint64_t s = lower.size();
        uniform_real_distribution<double> dist(0.0, 1.0);
        mt19937_64 mt;
        mcresult = monte_carlo_integration(s - 1, sampleNumber, trialNumber,
                                           ba, mt, dist, p);
        result.probability = mcresult.value;
        result.absoluteError = mcresult.error;
        result.relativeError = result.absoluteError / result.probability;
        result.success = true;
        result.upperBound = ba.getUpperBound();
        return result;
    }

}

#endif
