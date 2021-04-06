// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2020 Marios Papachristou

// Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_THIN_SAMPLES_HPP
#define DIAGNOSTICS_THIN_SAMPLES_HPP

template <typename NT, typename VT, typename MT>
MT thin_samples(MT const& samples, NT const& min_ess) {

    // Sample matrix is provided as d x n_samples
    unsigned int d = samples.rows();
    unsigned int N = samples.cols();
    unsigned int gap;
    unsigned int N_gap;

    // Thin samples are the initial samples which are N / min_ess apart
    gap = N / min_ess;
    N_gap = N - N % gap;

    MT thin_samples;
    thin_samples.resize(d, N_gap / gap);

    for (int i = 0; i < N_gap; i += gap) {
        thin_samples.col(i / gap) = samples.col(i);
    }

    return thin_samples;
}

#endif
