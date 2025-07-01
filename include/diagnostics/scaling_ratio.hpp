// VolEsti (volume computation and sampling library)

// Copyright (c) 

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_SCALING_RATIO_HPP
#define DIAGNOSTICS_SCALING_RATIO_HPP

template<typename Polytope, typename MT>
std::pair<std::vector<double>, typename Polytope::MT>scaling_ratio_boundary_test(const Polytope&           P,
                                                                                const MT&                 samples,
                                                                                const std::vector<int>&   facet_id,
                                                                                double                    tol       = 1e-10,
                                                                                double                    min_ratio = 0.01)
{
    using VT = typename Polytope::VT;
    const int dim = P.dimension();
    const int m = P.num_of_hyperplanes();
    const unsigned n_samp = static_cast<unsigned>(samples.cols());

    std::vector<double> scale(10); //we scale 10 times
    MT coverage(m, 10);

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
            double step = 0.1 * (k+1);
            double x = std::pow(step, 1.0 / dim);
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

    return { scale, coverage };
}

#endif // DIAGNOSTICS_SCALING_RATIO_HPP
