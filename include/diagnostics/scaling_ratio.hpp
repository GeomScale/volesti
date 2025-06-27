// VolEsti (volume computation and sampling library)

// Copyright (c) 

// Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef DIAGNOSTICS_SCALING_RATIO_HPP
#define DIAGNOSTICS_SCALING_RATIO_HPP

template< typename Polytope, typename MT, typename VT, typename StreamType >
void scaling_ratio_test(    const Polytope&   P,
                            const MT&         samples,
                            const VT&      facet_id,
                            double            tol,
                            StreamType&       out,
                            double            min_ratio = 0.01)
{
    const int true_dim = P.dimension();
    const int m = P.num_of_hyperplanes();
    const unsigned n_samp = static_cast<unsigned>(samples.cols());

    out << "\nScaling coverage by facet (skip < " << min_ratio << " ):\n";

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
        if (ratio < min_ratio) 
        {
            out << "facet " << f << " skipped (" << ratio << ")\n";
            continue;
        }


        // Finding the center
        Eigen::VectorXd p = Eigen::VectorXd::Zero(true_dim);
        for (int idx : S) p += samples.col(idx);
        p /= static_cast<double>(S.size());

        out << "Facet " << f << " (" << S.size() << " pts): ";

        // Looping over scale factor
        for (int k = 1; k <= 10; ++k) 
        {
            double step = 0.1 * k;
            double x = std::pow(step, 1.0 / true_dim);

            // Local copy of polytope for each scaling
            Polytope P_loc = P;
            
            //Shifting and scaling
            P_loc.shift(p);
            const Eigen::MatrixXd T = (1.0 / x) * Eigen::MatrixXd::Identity(true_dim, true_dim);
            P_loc.linear_transformIt(T);

            //Parameters of new polytope
            const auto& A_sh = P_loc.get_mat();
            const auto& b_sh = P_loc.get_vec();

            // Points still in facet
            unsigned survivors = 0;
            for (int idx : S) {

                const Eigen::VectorXd q_shift = samples.col(idx) - p;
                bool inside = true;

                for (int j = 0; j < A_sh.rows(); ++j) {
                    if (j == f) continue;
                    if (A_sh.row(j).dot(q_shift) - b_sh[j] > tol) {
                        inside = false; break;
                    }
                }
                if (inside) ++survivors;
            }

            double coverage = double(survivors) / double(S.size());
            out << std::pow(x, true_dim) << ':' << coverage;
            if (k < 10) out << ", ";

        }
        out << "\n";
    }
}


#endif 