#pragma once
#ifndef MCQMC_INTEGRATION_BOTEV_ALGORITHM_H
#define MCQMC_INTEGRATION_BOTEV_ALGORITHM_H
/**
 * @file BotevAlgorithm.h
 *
 * @brief A collection of functions to deal with the truncated
 * univariate and multivariate normal distributions.
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
 *
 * @note ported from
 * https://cran.r-project.org/package=TruncatedNormal
 * under License of GPL-3(https://cran.r-project.org/web/licenses/GPL-3)
 */

#include <vector>
#include <string>
#include <inttypes.h>

namespace MCQMCIntegration {
   /**
    * a class to compute Multivariate Normal Distribution.
    *
    * see:
    * @li https://cran.r-project.org/package=TruncatedNormal
    * @li Z. I. Botev (2015), The Normal Law Under Linear Restrictions:
    * Simulation and Estimation via Minimax Tilting, submitted to JRSS(B)
    * @li Z. I. Botev and P. L’Ecuyer (2015), Efficient Estimation and
    * Simulation of the Truncated Multivariate Student-t Distribution,
    * Proceedings of the 2015 Winter Simulation Conference, Huntington
    * Beach, CA, USA
    * @li Gibson G. J., Glasbey C. A., Elston D. A. (1994), Monte
    * Carlo evaluation of multivariate normal integrals and
    * sensitivity to variate ordering, In: Advances in Numerical
    * Methods and Applications, pages 120–126
    */
    class BotevAlgorithm {
    public:

        /**
         * a constructor of BotevAlgorithm from predefined parameter.
         *
         * @param[in] param_id predefined parameter ID.
         * @param[in] s_dimR dimension of R.
         */
        BotevAlgorithm(uint32_t param_id, uint32_t s_dimR);

        /**
         * a constructor of BotevAlgorithm.
         *
         * computes an estimator and a deterministic upper bound of the
         * probability P r(@c lower < X < @c upper), where X is a
         * zero-mean multivariate normal vector with covariance matrix Σ,
         * that is, X is drawn from N(0, Σ) infinite values for vectors @c
         * upper and @c lower are accepted;
         *
         * @param[in] lower lower truncation limit.
         * @param[in] upper upper truncation limit.
         * @param[in] sigma covariance matrix of N(0, Σ)
         */
        BotevAlgorithm(const std::vector<double>& lower,
                       const std::vector<double>& upper,
                       const std::vector< std::vector<double> >& sigma);

        /**
         * check and prepare.
         *
         * this member function should be called before Quasi
         * Monte-Carlo or Monte-Carlo integration.
         * @return true if parameter check and preparation success.
         */
        bool checkAndPrepare();

        /**
         * call back function for Quasi Monte-Carlo or Monte-Carlo integration.
         *
         * @param[in] x point in s -1 dimension, x should have length
         * s - 1, where s is length of @c lower and @c upper.
         * @return estimated value of a point @c x.
         */
        double operator()(const double x[]);

        /**
         * get theoretical upper bound on true Pr(@c lower < X < @c upper).
         * @return theoretical upper bound.
         */
        double getUpperBound() const;

        /**
         * destructor
         */
        ~BotevAlgorithm();

        /**
         * get number of predefined parameters.
         * @return number of predefined parameters.
         */
        static uint32_t getParameterSize();

        /**
         * get name of predefined parameter.
         * @param[in] index index of predefined parameter.
         * @return parameter name.
         */
        static const std::string getParameterName(uint32_t index);

        /**
         * get an explanation of predefined parameter.
         * @param[in] index index of predefined parameter.
         * @return an explanation of predefined parameter.
         */
        static const std::string getParameterConstruction(uint32_t index);
    private:
        // First of all, forbid copy and assign.
        BotevAlgorithm(const BotevAlgorithm& that);
        BotevAlgorithm& operator=(const BotevAlgorithm& that);
        class Impl;
        Impl *impl;
    };
}

#endif // MCQMC_INTEGRATION_BOTEV_ALGORITHM_H
