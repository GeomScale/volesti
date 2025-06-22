// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2018-2025 Apostolos Chalkis
// Copyright (c) 2025 Vladimir Necula

// Contributed and/or modified by Vladimir Necula, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "sampling/random_point_generators.hpp"
#include "random_walks/random_walks.hpp"
#include "convex_bodies/hpolytope.h"
#include "generators/h_polytopes_generator.h"
#include "generators/known_polytope_generators.h"
#include "diagnostics/effective_sample_size.hpp"
#include "diagnostics/multivariate_psrf.hpp"
#include "diagnostics/univariate_psrf.hpp"
#include "generators/order_polytope_generator.h"
#include "preprocess/barrier_center_ellipsoid.hpp"
#include "preprocess/svd_rounding.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

typedef HPolytope<Point> DenseHPOLYTOPE;
typedef HPolytope<Point, Eigen::SparseMatrix<NT, Eigen::RowMajor>> SparseHPOLYTOPE;
typedef BilliardWalk::template Walk<DenseHPOLYTOPE, RNGType> DenseBilliardWalkType;
typedef BilliardWalk::template Walk<SparseHPOLYTOPE, RNGType> SparseBilliardWalkType;

PushBackWalkPolicy push_back_policy;

const unsigned int FIXED_SEED = 42;

struct BenchmarkResults {
    NT ess_min;
    NT ess_avg;
    NT psrf_max;
    double time_walk;
    std::string walk_type;
    int dimension;
    int num_samples;
};

BenchmarkResults benchmark_dense_rounded_billiard_walk(DenseHPOLYTOPE& P, unsigned int num_samples, unsigned int walk_length) {
    std::cout << "Benchmarking Dense Rounded Billiard Walk..." << std::endl;

    using clock = std::chrono::high_resolution_clock;
    using seconds = std::chrono::duration<double>;

    RNGType rng(P.dimension());
    rng.set_seed(FIXED_SEED);

    auto [H, x_ac_vec, converged] = barrier_center_ellipsoid_linear_ineq<MT, EllipsoidType::LOG_BARRIER, NT>(P.get_mat(), P.get_vec());
    if (!converged) throw std::runtime_error("Failed to compute analytic center");
    Point x_ac(x_ac_vec);

    DenseHPOLYTOPE P_shifted = P;
    P_shifted.shift(-x_ac.getCoefficients());

    // Step 3: Compute the rounding transformation using Cholesky decompositio
    Eigen::LLT<MT> llt(H);
    if (llt.info() != Eigen::Success) throw std::runtime_error("LLT decomposition failed");
    
    MT L = llt.matrixL();
    MT L_inv = L.triangularView<Eigen::Lower>().solve(MT::Identity(P.dimension(), P.dimension()));
    
    MT A_rounded = P_shifted.get_mat() * L_inv;
    DenseHPOLYTOPE P_rounded(P.dimension(), A_rounded, P_shifted.get_vec());

    Point origin(P.dimension());
    origin.set_to_origin();
    if (!P_rounded.is_in(origin)) {
        origin = P_rounded.ComputeInnerBall().first;
    }

    auto t1 = clock::now();

    std::vector<Point> randPoints;
    typedef RandomPointGenerator<DenseBilliardWalkType> Generator;
    Generator::apply(P_rounded, origin, num_samples, walk_length, randPoints, push_back_policy, rng);

    auto t2 = clock::now();

    for (auto& p : randPoints) {
        VT y = p.getCoefficients();
        VT x_original = L_inv * y + x_ac.getCoefficients();
        p = Point(x_original);
    } 

    MT samples(P.dimension(), num_samples);
    for (size_t i = 0; i < randPoints.size(); ++i) {
        samples.col(i) = randPoints[i].getCoefficients();
    }

    NT psrf = multivariate_psrf<NT, VT, MT>(samples);
    unsigned int min_ess;
    VT ess_vector = effective_sample_size<NT, VT, MT>(samples, min_ess);

    BenchmarkResults results;
    results.ess_min = ess_vector.minCoeff();
    results.ess_avg = ess_vector.mean();
    results.psrf_max = psrf;
    results.time_walk = seconds(t2 - t1).count();
    results.walk_type = "Dense Rounded Billiard";
    results.dimension = P.dimension();
    results.num_samples = num_samples;

    return results;
}


BenchmarkResults benchmark_sparse_billiard_walk(SparseHPOLYTOPE& P, unsigned int num_samples, unsigned int walk_length) {
    std::cout << "Benchmarking Sparse Billiard Walk..." << std::endl;
    
    using clock = std::chrono::high_resolution_clock;
    using seconds = std::chrono::duration<double>;
    
    RNGType rng(P.dimension());
    rng.set_seed(FIXED_SEED);
    
    auto A_matrix = P.get_mat();
    auto b_vector = P.get_vec();

    MT Hessian;
    VT x_ac_vec;
    bool converged;

    MT A_dense = MT(A_matrix);
    VT b_dense = b_vector;

    auto result = barrier_center_ellipsoid_linear_ineq<MT, EllipsoidType::LOG_BARRIER, NT>(A_dense, b_dense);

    Hessian = std::get<0>(result);
    x_ac_vec = std::get<1>(result);
    converged = std::get<2>(result); 
    
    Point x_ac(x_ac_vec);
    
    SparseHPOLYTOPE P_shifted = P;
    P_shifted.shift(-x_ac.getCoefficients());
    
    Eigen::SparseMatrix<NT, Eigen::ColMajor> H_sparse;
    if constexpr (std::is_same_v<decltype(Hessian), Eigen::SparseMatrix<NT, Eigen::ColMajor>>) {
        H_sparse = Hessian;
    } else {
        H_sparse = Hessian.sparseView();
    }
    
    Point origin(P.dimension());
    origin.set_to_origin();

    if ( !P_shifted.is_in(origin) ) {
        origin = P_shifted.ComputeInnerBall().first;
        std::cout << "Computed inner ball" << std::endl;
    } else {
        std::cout << "WARNING: Origin is not in the polytope" << std::endl;
        return BenchmarkResults();
    }
    
    typedef SparseBilliardWalk::template Walk<SparseHPOLYTOPE, RNGType> SparseBilliardWalkType;
    SparseBilliardWalk::parameters parms(walk_length, true);
    SparseBilliardWalk::Walk<SparseHPOLYTOPE, RNGType> walk(P_shifted, origin, rng, parms, H_sparse);

    auto t1 = clock::now();

    std::vector<Point> randPoints;
    randPoints.reserve(num_samples);

    Point current_point = origin;
    unsigned int successful_samples = 0;

    for (unsigned int i = 0; i < num_samples; ++i) {
        try {
            walk.apply(P_shifted, current_point, walk_length, rng);
            randPoints.push_back(Point(current_point.getCoefficients() + x_ac.getCoefficients()));
            successful_samples++;
        } catch (const std::exception& e) {
            std::cerr << "Error during walk step " << i << ": " << e.what() << std::endl;
            continue;
        }
    } 
    
    auto t2 = clock::now();
    
    if (successful_samples == 0) {
        throw std::runtime_error("No successful samples were generated");
    }
    
    MT samples(P.dimension(), successful_samples);
    for (size_t i = 0; i < randPoints.size(); ++i) {
        samples.col(i) = randPoints[i].getCoefficients();
    }
    
    NT psrf = multivariate_psrf<NT, VT, MT>(samples);
    unsigned int min_ess;
    VT ess_vector = effective_sample_size<NT, VT, MT>(samples, min_ess);
    
    BenchmarkResults results;
    results.ess_min = ess_vector.minCoeff();
    results.ess_avg = ess_vector.mean();
    results.psrf_max = psrf;
    results.time_walk = seconds(t2 - t1).count();
    results.walk_type = "Sparse Billiard";
    results.dimension = P.dimension();
    results.num_samples = successful_samples;
    
    return results;
}


void print_results(const std::vector<BenchmarkResults>& results) {
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "BENCHMARK RESULTS" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    std::cout << std::left << std::setw(15) << "Walk Type"
              << std::setw(8) << "Dim " 
              << std::setw(10) << "Samples "
              << std::setw(12) << "Time Walk(s) "
              << std::setw(12) << "Min ESS "
              << std::setw(12) << "Avg ESS "
              << std::setw(12) << "Max PSRF " << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    for (const auto& result : results) {
        std::cout << std::left << std::setw(15) << result.walk_type
                  << std::setw(8) << result.dimension
                  << std::setw(10) << result.num_samples
                  << std::setw(12) << std::fixed << std::setprecision(3) << result.time_walk
                  << std::setw(12) << std::fixed << std::setprecision(1) << result.ess_min
                  << std::setw(12) << std::fixed << std::setprecision(1) << result.ess_avg
                  << std::setw(12) << std::fixed << std::setprecision(3) << result.psrf_max << std::endl;
    }
    std::cout << std::string(80, '=') << std::endl; 
}

void run_benchmark_case(std::vector<BenchmarkResults>& all_results, 
                       unsigned int dim, 
                       unsigned int num_relations,
                       unsigned int num_samples,
                       const std::string& test_name) {
    std::cout << "\n=== " << test_name << " ===" << std::endl;
    
    DenseHPOLYTOPE P_dense = random_orderpoly<DenseHPOLYTOPE, NT>(dim, num_relations, FIXED_SEED);
    P_dense.ComputeInnerBall();
    
    SparseHPOLYTOPE P_sparse(P_dense.dimension(), P_dense.get_mat().sparseView(), P_dense.get_vec());
    P_sparse.ComputeInnerBall();
    
    NT walk_length = NT(1) * dim; 
    
    auto dense_rounded_result = benchmark_dense_rounded_billiard_walk(P_dense, num_samples, walk_length);
    auto sparse_result = benchmark_sparse_billiard_walk(P_sparse, num_samples, walk_length);
    
    all_results.push_back(dense_rounded_result);
    all_results.push_back(sparse_result);
}


void run_comprehensive_benchmark() {
    std::vector<BenchmarkResults> all_results;
    unsigned int num_samples = 5000;
    
    // Test 1: Small order polytope (10D) with sparse relations
    run_benchmark_case(all_results, 10, 25, num_samples, "Test 1: 10D Order Polytope (Sparse)");
    
    // Test 2: Medium order polytope (15D) with medium density
    run_benchmark_case(all_results, 15, 45, num_samples, "Test 2: 15D Order Polytope (Medium)");
    
    // Test 3: Large order polytope (20D) with dense relations
    run_benchmark_case(all_results, 20, 80, num_samples, "Test 3: 20D Order Polytope (Dense)");
    
    // Test 4: High-dimensional polytope (30D) with sparse relations
    run_benchmark_case(all_results, 30, 100, num_samples, "Test 4: 30D Order Polytope (Sparse)");
    
    // Test 5: Very high-dimensional polytope (40D) with medium density
    run_benchmark_case(all_results, 40, 150, num_samples, "Test 5: 40D Order Polytope (Medium)");
    
    // Test 6: Ultra high-dimensional polytope (50D) with sparse relations
    run_benchmark_case(all_results, 50, 200, num_samples, "Test 6: 50D Order Polytope (Sparse)");
    
    print_results(all_results);
}

int main() {
    std::cout << "Sparse Billiard Walk Benchmark" << std::endl;
    std::cout << "===============================" << std::endl;
    std::cout << "This benchmark compares regular vs sparse billiard walk performance" << std::endl;
    std::cout << "using Effective Sample Size (ESS) and PSRF metrics." << std::endl;
    try {
        run_comprehensive_benchmark();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (const char* msg) {
        std::cerr << "Error: " << msg << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
    return 0;
} 