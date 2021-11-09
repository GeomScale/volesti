// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_PSRF_UPDATER_HPP
#define DIAGNOSTICS_PSRF_UPDATER_HPP

//#include "ess_updater_autocovariance.hpp"


/**
   This is a class that updates the effective sample size (ess) of a sample given a new chain
   using Welford's algorithm to update the average values and the variance estimates where needed.
   The chains has to be of the same length. The ess estimation exploits Geyer's stable estimator
   for the autocovariance and the Geyer's conversion to a monotone sequence, given in,

   Charles J. Geyer, Practical Markov Chain Monte Carlo, Statistical Science 1992.

 * @tparam NT number type
 * @tparam VT vector type
 * @tparam MT matrix type
*/
template <typename NT, typename VT, typename MT>
class PSRFestimator {

private:
   unsigned int         middle_indx, window_len, d, middle_indx_prev, num_samples; 
   VT                   mean1, mean2, mean00, sigma1, sigma2, R, W, B, sigma, sum_sigma, M2, M1, delta1, delta2, sum2, sum_sq2, delta0, temp1, temp2;
   //NT                   oldM, rho_hat_odd, rho_hat_even, mean_var, new_elem;
   MT                   acov_s_mean, rho_hat_s;
   
public:
    PSRFestimator() {}

    PSRFestimator(unsigned int const& _window_len, unsigned int const& _dim) 
    {
        window_len = _window_len;
        d = _dim;
        num_samples = 0;
        middle_indx = 0;
        middle_indx_prev = 0;

        mean1.setZero(d);
        mean2.setZero(d);
        mean00.setZero(d);
        sigma1.setZero(d);
        sigma2.setZero(d);
        R.setZero(d);
        W.setZero(d);
        B.setZero(d); 
        sigma.setZero(d);
        sum_sigma.setZero(d);
        M1.setZero(d);
        delta1.setZero(d);
        M2.setZero(d);
        delta2.setZero(d);
        delta0.setZero(d);
        sum2.setZero(d);
        sum_sq2.setZero(d);
        temp1.setZero(d);
        temp2.setZero(d);
        //draws.setZero(num_draws);
        //acov_s_mean.setZero(num_draws-3, d);
        //rho_hat_s.setZero(num_draws, d);
   }

    void update_estimator(MT const& samples) 
    {
        num_samples += window_len;
        middle_indx += window_len / 2;
        int counter = 0;
   
        for (int i = middle_indx_prev; i < middle_indx; i++)
        {
            temp1 = samples.col(i);
            temp2 = samples.col(num_samples - window_len + counter);

            delta1 = temp1 - mean1;
            mean1 += delta1 / (i + 1);
            M1 += delta1.cwiseProduct(temp1 - mean1);
            sigma1 = M1 / (i + 1);

            delta2 = temp2 - mean2;
            mean2 += delta2 / (i + 1);
            M2 += delta2.cwiseProduct(temp2 - mean2);
            sigma2 = M2 / (i + 1);
            sum2 += temp2;
            sum_sq2 += temp2.cwiseProduct(temp2);

            delta0 = temp2 - mean00;
            mean00 += delta0 / (counter + num_samples - window_len + 1);

            counter++;
        }

        for (int i = 0; i < (window_len / 2); i++)
        {
            temp2 = samples.col(num_samples - window_len/2 + i);

            delta2 = temp2 - mean2;
            mean2 += delta2 / (middle_indx + i + 1);
            M2 += delta2.cwiseProduct(temp2 - mean2);
            sigma2 = M2 / (middle_indx + i + 1);

            sum2 += temp2;
            sum_sq2 += temp2.cwiseProduct(temp2);

            delta0 = temp2 - mean00;
            mean00 += delta0 / (i + num_samples - window_len/2 + 1);
        }

        unsigned int len = middle_indx + window_len / 2;
        for (int i = middle_indx_prev; i < middle_indx; i++)
        {
            temp2 = samples.col(i);

            mean2 = sum2 / (len - 1) - temp2 / (len - 1); 
            sum2 -= temp2;
            sum_sq2 -= temp2.cwiseProduct(temp2);

            len--;
        }

        sigma2 = (sum_sq2 - NT(2) * (sum2.cwiseProduct(mean2)) + NT(len) * mean2.cwiseProduct(mean2)) * (NT(1) / (NT(len) - NT(1)) );

        middle_indx_prev = middle_indx;
    }


    void estimate_psrf()
    {
        W = (sigma1 + sigma2) / NT(2);
        temp1 = mean1 - mean00;
        temp2 = mean2 - mean00;

        B = temp1.cwiseProduct(temp1) + temp2.cwiseProduct(temp2);
        sigma = ((NT(middle_indx) - NT(1)) / NT(middle_indx)) * W + B;

        NT* sigma_data = sigma.data();
        NT* W_data = W.data();
        NT* R_data = R.data();

        for (int i=0; i<d; i++)
        {
            *R_data = std::sqrt((*sigma_data) / (*W_data));

            R_data++;
            sigma_data++;
            W_data++;
        }

        //R = std::sqrt(sigma / W);
    }


    VT get_psrf() 
    {
      return R;
    }

};


#endif

