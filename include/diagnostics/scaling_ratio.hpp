// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025-2025 Iva Janković

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_SCALING_RATIO_HPP
#define DIAGNOSTICS_SCALING_RATIO_HPP

template<typename Polytope>
std::tuple<typename Polytope::VT, typename Polytope::MT,typename Polytope::VT, typename Polytope::VT>
scaling_ratio_boundary_test(const Polytope&  P,
                            const typename Polytope::MT& samples,
                            const typename Polytope::NT tol = 1e-10,
                            const typename Polytope::NT min_ratio = 0.01) 
{
    using VT = typename Polytope::VT;
    using MT = typename Polytope::MT;
    using NT = typename Polytope::NT;
    
    const int dim = P.dimension();
    const int m = P.num_of_hyperplanes();
    const unsigned n_samp = static_cast<unsigned>(samples.cols());

    VT scale(10); //we scale 10 times
    MT coverage(m, 10);
    std::vector<int> facet_id(n_samp, -1);
    const auto A_full = P.get_mat();
    const auto b_full = P.get_vec();

    for (int i = 0; i < n_samp; ++i) {
        auto Aq = A_full * samples.col(i);
        
        for (size_t k = 0; k < m; ++k) {
            if (std::abs(Aq[k] - b_full[k]) < tol) 
            {
                facet_id[i] = static_cast<int>(k);
                break;
            }
        }
    }

    for (int f = 0; f < m; ++f) 
    {
        // Samples S on facet f
        std::vector<int> S;
        S.reserve(n_samp / m);
        for (unsigned i = 0; i < n_samp; ++i) {
            if (facet_id[i] == f)
                S.push_back(static_cast<int>(i));
        }

        const double ratio = static_cast<double>(S.size()) / n_samp;
        if (ratio < min_ratio)  continue;

        // Finding the center
        VT p = VT::Zero(dim);
        for (int idx : S) p += samples.col(idx);
        p /= static_cast<double>(S.size());

        // Looping over scale factor
        for (int k = 0; k < 10; ++k) 
        {
            NT step = 0.1 * (k+1);
            NT x = std::pow(step, 1.0 / dim);
            // Local copy of polytope for each scaling
            Polytope P_loc = P;
            
            //Shifting and scaling
            P_loc.shift(p);
            MT T = (1.0 / x) * MT::Identity(dim, dim);
            P_loc.linear_transformIt(T);

            //Parameters of new polytope
            const auto& A_sh = P_loc.get_mat();
            const auto& b_sh = P_loc.get_vec();

            // Points still in facet
            unsigned survivors = 0;
            for (int idx : S) {

                const VT q_shift = samples.col(idx) - p;
                bool inside = true;

                for (int j = 0; j < A_sh.rows(); ++j) {
                    if (j == f) continue;
                    if (A_sh.row(j).dot(q_shift) - b_sh[j] > tol) {
                        inside = false; break;
                    }
                }
                if (inside) ++survivors;
            }

            coverage(f, k) = double(survivors) / double(S.size());
            scale[k]=std::pow(x, dim);
        }
    }

    int K = scale.size();           // number of scaling factors
    VT max_dev(m), avg_dev(m);      // Vectors of average and maximum deviations

    for (int f = 0; f < m; ++f) {
        double sumd = 0.0;
        double maxd = 0.0;
        for (int k = 0; k < K; ++k) {
            double d = std::abs(coverage(f, k) - scale[k])* 100.0; // in percentage
            sumd += d;
            if (d > maxd) maxd = d;
        }
        avg_dev[f] = sumd/K;
        max_dev[f] = maxd;      
    }

    return {scale, coverage, max_dev, avg_dev};
}

#endif // DIAGNOSTICS_SCALING_RATIO_HPP
