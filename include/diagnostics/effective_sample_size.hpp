// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#include <unsupported/Eigen/FFT>

#ifndef EFFECTIVE_SAMPLE_SIZE_HPP
#define EFFECTIVE_SAMPLE_SIZE_HPP

template <typename NT>
void cummulative_minimum(std::vector<NT> &v) {
    unsigned int N = v.size();
    for (unsigned int i = 1; i < N; i++) {
        if (v[i] > v[i - 1]) v[i] = v[i - 1];
    }
}

template <typename NT, typename VT, typename MT>
VT effective_sample_size(MT const& samples, unsigned int &min_ess) {
    typedef Eigen::FFT<NT> EigenFFT;
    typedef std::complex<NT> CNT;
    EigenFFT fft;

    // Sample matrix is provided as d x n_samples
    unsigned int d = samples.rows();
    unsigned int N = samples.cols();
    unsigned int N_even = N - N % 2;

    // Calculate effective sample size per dimension
    VT ess;
    ess.resize(d);

    // Autocorrelation vector
    std::vector<NT> min_auto_correlation(N_even / 2, NT(0));


    // Z-normalized samples
    std::vector<NT> normalized_sample_row(N);

    // FFT vector
    std::vector<CNT> fft_vec(N);
    std::vector<CNT> fft_inv_vec(N);
    std::vector<CNT> psd(N);

    // Helper variables
    NT row_mean;
    NT variance;
    NT gap;

    for (unsigned int i = 0; i < d; i++) {

        // Z-normalization
        row_mean = samples.row(i).mean();

        for (int j = 0; j < N; j++) {
            normalized_sample_row[j] = samples(i, j) - row_mean;
        }

        variance = NT(0);

        for (int j = 0; j < N; j++) {
            variance += pow(normalized_sample_row[j], 2);
        }

        variance *= (1.0 / N);

        for (int j = 0; j < N; j++) {
            normalized_sample_row[j] /= NT(1e-16 + sqrt(variance));
        }

        // Perform FFT on 2N points
        fft.fwd(fft_vec, normalized_sample_row);

        // Calculate PSD which is the norm squared of the FFT of the sequence
        for (int j = 0; j < N; j++) {
            psd[j].real(std::norm(fft_vec[j]));
            psd[j].imag(NT(0));
        }

        // Invert fft to get autocorrelation function
        fft.inv(fft_inv_vec, psd);

        // Calculate minimum autocorrelation
        for (int j = 0; j < N_even; j += 2) {
            min_auto_correlation[j / 2] =
                (1.0 / N) * (std::real(fft_inv_vec[j + 1]) +
                std::real(fft_inv_vec[j]));
        }

        gap = NT(0);
        cummulative_minimum(min_auto_correlation);

        for (int j = 0; j < N_even / 2; j++) {
            if (min_auto_correlation[j] > 0) gap += min_auto_correlation[j];
        }

        gap = 2 * gap - NT(1);
        if (gap < NT(1)) gap = NT(1);

        ess(i) = (1.0 * N) / gap;

        // Store minimum effective sample size as integer (in general ess is not int)
        // for the thinning process
        if (i == 0) min_ess = (unsigned int)ess(i);
        else if ((unsigned int)ess(i) < min_ess) min_ess = (unsigned int)ess(i);

    }

    return ess;
}

#endif
