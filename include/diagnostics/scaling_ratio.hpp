// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_SCALING_RATIO_HPP
#define DIAGNOSTICS_SCALING_RATIO_HPP

template<typename Polytope>
std::tuple<typename Polytope::VT,typename Polytope::MT,typename Polytope::VT,typename Polytope::VT,int,double>
scaling_ratio_boundary_test(const Polytope& P,const typename Polytope::MT& samples,const typename Polytope::NT tol = 1e-10,const typename Polytope::NT min_ratio = 0.01)
{
    using VT = typename Polytope::VT;
    using MT = typename Polytope::MT;
    using NT = typename Polytope::NT;

    /*
    DISCLAIMER (Keep in mind when interpreting results):
    This diagnostic assumes that P is a full-dimensional H-polytope.
    The facet identification and survivor logic rely on the assumption
    that each facet corresponds to a unique supporting hyperplane and that
    removing one constraint locally produces a (dim-1)-dimensional facet.
    For degenerate polytopes the results may be undefined.
    */

    const int dim = P.dimension();
    const int m   = P.num_of_hyperplanes();
    const int K   = 10;
    const unsigned n_samp = static_cast<unsigned>(samples.cols());

    // Always define scale
    VT scale(K);
    for (int k = 0; k < K; ++k)
        scale[k] = NT(0.1) * NT(k + 1);

    // Pre-initialize outputs (deterministic, no UB)
    MT coverage(m, K);
    coverage.setConstant(std::numeric_limits<double>::quiet_NaN());

    VT max_dev(m), avg_dev(m);
    max_dev.setConstant(std::numeric_limits<double>::quiet_NaN());
    avg_dev.setConstant(std::numeric_limits<double>::quiet_NaN());

    int zero_facets_count = 0;
    double zero_facets_percent = 0.0;

    // Guard: no samples
    if (n_samp == 0) {
        scale.setZero();
        zero_facets_count = m;
        zero_facets_percent = (m > 0) ? 100.0 : 0.0;
        return {scale, coverage, max_dev, avg_dev,
                zero_facets_count, zero_facets_percent};
    }

    // Guard: low dimension
    if (dim <= 1) {
        return {scale, coverage, max_dev, avg_dev,
                zero_facets_count, zero_facets_percent};
    }

    std::vector<int> facet_id(n_samp, -1);

    const auto A_full = P.get_mat();
    const auto b_full = P.get_vec();

    // Assign each sample to first matching facet
    for (unsigned i = 0; i < n_samp; ++i) {
        auto Aq = A_full * samples.col(static_cast<int>(i));
        for (int k = 0; k < m; ++k) {
            if (std::abs(Aq[k] - b_full[k]) < tol) {
                facet_id[i] = k;
                break;
            }
        }
    }

    // Per-facet loop
    for (int f = 0; f < m; ++f) {

        std::vector<int> S;
        S.reserve(std::max<unsigned>(1u,n_samp / static_cast<unsigned>(std::max(1, m))));

        for (unsigned i = 0; i < n_samp; ++i)
            if (facet_id[i] == f)
                S.push_back(static_cast<int>(i));

        if (S.empty())
            zero_facets_count++;

        const double ratio =static_cast<double>(S.size()) /static_cast<double>(n_samp);

        if (S.empty() || ratio < min_ratio)
            continue;  // leave NaNs

        // Compute facet center
        VT p = VT::Zero(dim);
        for (int idx : S)
            p += samples.col(idx);
        p /= static_cast<double>(S.size());

        for (int k = 0; k < K; ++k) {

            const NT step = scale[k];
            const NT x = std::pow(step, NT(1) / NT(dim));

            Polytope P_loc = P;
            P_loc.shift(p);

            MT T = (NT(1) / x) * MT::Identity(dim, dim);
            P_loc.linear_transformIt(T);

            const auto& A_sh = P_loc.get_mat();
            const auto& b_sh = P_loc.get_vec();

            // Find corresponding facet plane
            int f_loc = -1;
            NT best = std::numeric_limits<NT>::infinity();
            for (int j = 0; j < A_sh.rows(); ++j) {
                NT v = std::abs(b_sh[j]);
                if (v < best) {
                    best = v;
                    f_loc = j;
                }
            }

            const bool have_f_loc =(f_loc >= 0 && best < NT(100) * tol);

            unsigned survivors = 0;

            for (int idx : S) {

                const VT q_shift =samples.col(idx) - p;

                bool inside = true;

                for (int j = 0; j < A_sh.rows(); ++j) {
                    if (have_f_loc && j == f_loc)
                        continue;

                    if (A_sh.row(j).dot(q_shift)
                        - b_sh[j] > tol) {
                        inside = false;
                        break;
                    }
                }

                if (inside)
                    ++survivors;
            }

            coverage(f, k) =static_cast<double>(survivors) /static_cast<double>(S.size());
        }
    }

    zero_facets_percent =(m > 0)? (100.0 * static_cast<double>(zero_facets_count)/ static_cast<double>(m)): 0.0;

    // Deviations
    for (int f = 0; f < m; ++f) {
        double sumd = 0.0;
        double maxd = 0.0;
        int cnt = 0;
        for (int k = 0; k < K; ++k) {
            const double c = coverage(f, k);
            if (!std::isfinite(c))
                continue;

            const double dperc =std::abs(c -static_cast<double>(scale[k])) * 100.0;

            sumd += dperc;
            if (dperc > maxd)
                maxd = dperc;
            ++cnt;
        }
        if (cnt > 0) {
            avg_dev[f] = sumd /
                static_cast<double>(cnt);
            max_dev[f] = maxd;
        }
    }

    return {scale, coverage, max_dev, avg_dev,zero_facets_count, zero_facets_percent};

}

#endif // DIAGNOSTICS_SCALING_RATIO_HPP