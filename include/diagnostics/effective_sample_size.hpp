// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements the Rubin & Gelman diagnostic.
    It is based on "Inference from iterative simulation using multiple sequences, 1992" by D. B. Rubin and A. Gelman

    For each coordinate the sample is splitted into two parts. 
    Then the psrf of D.B. Rubin and A. Gelman is computed for each coordinate

    Inputs: samples, a matrix that contains sample points column-wise

    Output: The value of PSRF of D.B. Rubin and A. Gelman for each coordinate
*/

#ifndef EFFECTIVE_SAMPLE_SIZE_HPP
#define EFFECTIVE_SAMPLE_SIZE_HPP

template <typename NT, typename VT, typename MT>
VT eff_univariate(MT samples)
{
    //MT runs = samples.transpose();
    unsigned int N = samples.cols(), d = samples.rows();
    VT coord_samples(N), results(d);
    NT st_dev, numenator, sum, rhat_prev, rhat_curr, curr_sum, prev_sum, sum_rhats;

    //VT mean_mat = samples.rowwise().mean();

    samples.noalias() = samples.colwise() - samples.rowwise().mean();

    for (int i = 0; i < d; i++)
    {
        coord_samples = samples.row(i);
        curr_sum = 0.0;
        prev_sum = 0.0;
        sum_rhats = 0.0;
        //rhats.setZero(N-1);

        for (int k = 1; k < N; k++) 
        {
            sum = 0.0;
            for (int j = 0; j < N; j++)
            {
                sum += coord_samples.coeff(j) * coord_samples.coeff(j);
            }
            
            st_dev = (1.0 / (NT(N) - 1.0)) * sum;
            sum = 0.0;
            for (int j = 0; j < (N-k); j++)
            {
                sum += coord_samples.coeff(j) * coord_samples.coeff(k + j);
            }
            
            numenator = (1.0 / NT(N)) * sum;
            rhat_curr = numenator / st_dev;
            
            if (k%2 == 1 && k > 1) {
                curr_sum = rhat_prev + rhat_curr;
                if (curr_sum > 0.0) {
                    if (curr_sum > prev_sum) {
                        sum_rhats += (curr_sum);
                    } else {
                        //std::cout<<"[NO monotonic] k = "<<k<<", rhat_prev + rhat_curr = "<<rhat_prev + rhat_curr<<", rhat = "<<rhat<<std::endl;
                        break;
                    }
                } else {
                    //std::cout<<"[NO positive] k = "<<k<<", rhat_prev + rhat_curr = "<<rhat_prev + rhat_curr<<", rhat = "<<rhat<<std::endl;
                    break;
                }
                prev_sum = curr_sum;
            } else if (k == 1) {
                sum_rhats += rhat_curr;
            }
            rhat_prev = rhat_curr;
        }
        results(i) = int(NT(N) / (1.0 + 2.0 * sum_rhats));
    }
    return results;
}

#endif

