// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements the interval diagnostic of Brooks & Gelman.
    It is based on  "General Methods for Monitoring Convergence of Iterative Simulations, 1998" by S. Brooks and A. Gelman

    For each coordinate the sample is splitted into two parts.
    Then the psrf of S. Brooks and A. Gelman is computed for each coordinate.

    Inputs: samples, a matrix that contains sample points column-wise

    Output: The value of interval PSRF of S. Brooks and A. Gelman for each coordinate
*/

#ifndef DIAGNOSTICS_INTERVAL_PSRF_HPP
#define DIAGNOSTICS_INTERVAL_PSRF_HPP

template <typename VT, typename NT, typename MT>
VT interval_psrf(MT const& samples, NT alpha = 0.05)
{
    MT runs = samples.transpose();
    unsigned int N = samples.cols(), d = samples.rows();
    unsigned int N1 = N / 2;
    unsigned int N2 = N - N1;
    VT sorted_samples(N), marginal_samples(N), sorted_subsamples1(N1), sorted_subsamples2(N2), results(d);
    std::vector<NT> temp_col(N);

    for (int i = 0; i < d; i++)
    {
        sorted_samples = runs.col(i);
        marginal_samples = runs.col(i);

        temp_col.resize(N);
        temp_col = std::vector<NT>(&sorted_samples[0], sorted_samples.data() + sorted_samples.cols() *
                                   sorted_samples.rows());
        std::sort(temp_col.begin(), temp_col.end());
        sorted_samples = Eigen::Map<VT>(&temp_col[0], temp_col.size());

        int n1 = N * (alpha / NT(2)), n2 = N - N * (alpha / NT(2));

        NT len_total_sequence_interval = sorted_samples(n2) - sorted_samples(n1);

        sorted_subsamples1 = marginal_samples.block(0,0,N1,1);
        temp_col.resize(N1);
        temp_col = std::vector<NT>(&sorted_subsamples1[0], sorted_subsamples1.data() +
                                   sorted_subsamples1.cols() * sorted_subsamples1.rows());
        std::sort(temp_col.begin(), temp_col.end());
        sorted_subsamples1 = Eigen::Map<VT>(&temp_col[0], temp_col.size());

        sorted_subsamples2 = marginal_samples.block(N1,0,N2,1);
        temp_col.resize(N2);
        temp_col = std::vector<NT>(&sorted_subsamples2[0], sorted_subsamples2.data() +
                                   sorted_subsamples2.cols() * sorted_subsamples2.rows());
        std::sort(temp_col.begin(), temp_col.end());
        sorted_subsamples2 = Eigen::Map<VT>(&temp_col[0], temp_col.size());

        n1 = N1 * (alpha / NT(2)), n2 = N1 - N1 * (alpha / NT(2));
        NT len_sequence_interval1 = sorted_subsamples1(n2) - sorted_subsamples1(n1);

        n1 = N2 * (alpha / NT(2)), n2 = N2 - N2 * (alpha / NT(2));
        NT len_sequence_interval2 = sorted_subsamples2(n2) - sorted_subsamples2(n1);

        NT R = (len_total_sequence_interval) /
                ((len_sequence_interval1 + len_sequence_interval2) / NT(2));

        results(i) = std::abs(1.0 - R) + NT(1);
     }
     return results;
}


#endif
