// VolEsti (volume computation and sampling library)

//KS Test (radial distribution)

#ifndef DIAGNOSTICS_KS_TEST_HPP
#define DIAGNOSTICS_KS_TEST_HPP

#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// Kolmogorov distribution: P(K > z)
inline double kolmogorov_prob(double z) {
    if (z <= 0.0) return 1.0; 
    double sum = 0.0;
    for (int k = 1; k <= 200; ++k) {
        double term = std::exp(-2.0 * k * k * z * z);
        sum += (k % 2 ? term : -term);
        if (term < 1e-15) break; 
    }
    return std::max(0.0, std::min(1.0, 2.0 * sum));
}

// Global uniformity test
template<typename Polytope>
std::tuple<double, double, std::vector<double>, std::vector<double>>
global_scaling_test(const Polytope& P,
                    const typename Polytope::MT& samples,
                    int thinning_factor = 1) 
{
    using VT = typename Polytope::VT;

    const int dim    = P.dimension();
    const int n_total = static_cast<int>(samples.cols());

    // Setup center and constraints
    VT center = samples.rowwise().mean();
    const auto A = P.get_mat();   
    const auto b = P.get_vec();   
    VT b_shifted = b - A * center;

    if (b_shifted.minCoeff() < 1e-12) {
        std::cerr << "[GlobalKS] Warning: Empirical center is on boundary/outside.\n";
    }
    std::vector<double> rvals_all;
    rvals_all.reserve(n_total);

    for (int i = 0; i < n_total; ++i) {
        VT q = samples.col(i) - center;
        VT u = A * q;
        
        double r_max = 0.0;
        for (int k = 0; k < u.size(); ++k) {
            if (u[k] > 0.0) {
                double denom = b_shifted[k];
                if (denom > 1e-14) {
                    double t = u[k] / denom;
                    if (t > r_max) r_max = t;
                } else {
                    r_max = 1.0; 
                }
            }
        }
        if (r_max > 1.0) r_max = 1.0;
        if (r_max < 0.0) r_max = 0.0;
        rvals_all.push_back(r_max);
    }

    // We only perform the KS test on the subset to ensure independence.
    std::vector<double> uvals;
    if (thinning_factor < 1) thinning_factor = 1;
    
    // Reserve roughly N / factor
    uvals.reserve(n_total / thinning_factor + 1);

    for (int i = 0; i < n_total; i += thinning_factor) {
        double r = rvals_all[i];
        // Transform r -> u = r^d (Probability Integral Transform)
        uvals.push_back(std::pow(r, dim));
    }

    // Compute KS Statistic 
    std::sort(uvals.begin(), uvals.end());
    int N_test = static_cast<int>(uvals.size());

    if (N_test < 5) return {0.0, 1.0, {}, {}}; // Too few samples

    double ks_stat = 0.0;
    for (int i = 0; i < N_test; ++i) {
        double F_emp_lo = static_cast<double>(i) / N_test;
        double F_emp_hi = static_cast<double>(i + 1) / N_test;
        double F_theo   = uvals[i];

        double diff = std::max(std::abs(F_emp_lo - F_theo), 
                               std::abs(F_emp_hi - F_theo));
        if (diff > ks_stat) ks_stat = diff;
    }

    // Standard P-value 
    double sqrt_n = std::sqrt(static_cast<double>(N_test));
    double lambda = (sqrt_n + 0.12 + 0.11 / sqrt_n) * ks_stat;
    double p_value = kolmogorov_prob(lambda);

    // Shell diagnostics
    std::vector<double> exp_coverage(10), obs_coverage(10, 0.0);
    std::vector<double> r_thresholds(10);
    std::vector<int>    shell_counts(10, 0);

    for (int k = 0; k < 10; ++k) {
        exp_coverage[k] = 0.1 * (k + 1);
        r_thresholds[k] = std::pow(exp_coverage[k], 1.0 / dim);
    }

    for (double r : rvals_all) {
        for (int k = 0; k < 10; ++k) {
            if (r <= r_thresholds[k]) {
                shell_counts[k]++;
                break;
            }
        }
    }

    int cumulative = 0;
    for (int k = 0; k < 10; ++k) {
        cumulative += shell_counts[k];
        obs_coverage[k] = static_cast<double>(cumulative) / n_total;
    }

    return {ks_stat, p_value, obs_coverage, exp_coverage};
}

#endif
