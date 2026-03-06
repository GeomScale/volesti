#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <string>

#include <boost/random.hpp>
#include "Eigen/Eigen"

// These help with generating polytopes
#include "hpolytope.h"
#include "known_polytope_generators.h"
#include "custom_generators.h"
#include "order_polytope_generator.h"

#include "cartesian_geom/cartesian_kernel.h"
#include "sampling/random_point_generators.hpp"
#include "random_walks/random_walks.hpp"
#include "random_walks/sparse_uniform_billiard_walk.hpp"
#include "preprocess/max_inscribed_ellipsoid.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "convex_bodies/ellipsoid.h"
#include "convex_bodies/hpolytope.h"

// These are related to PSRF and ESS and K-S test
#include "sampling/sample_correlation_matrices.hpp"
#include "matrix_operations/EigenvaluesProblems.h"
#include "diagnostics/effective_sample_size.hpp"
#include "diagnostics/univariate_psrf.hpp"
#include "diagnostics/scaling_ratio.hpp"
#include "diagnostics/KS_test.hpp"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef HPolytope <Point> HPOLYTOPE;

//Usefull struct to rotate polytope
template <typename HPOLYTOPE>
HPOLYTOPE rotate_all_dims(const HPOLYTOPE& P, typename HPOLYTOPE::NT angle)
{
    using NT = typename HPOLYTOPE::NT;
    int dim = P.dimension();

    // Build global rotation matrix
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> R =
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Identity(dim, dim);

    NT c = std::cos(angle);
    NT s = std::sin(angle);

    // Apply rotation in each adjacent coordinate plane
    for (int k = 0; k < dim - 1; ++k) {
        Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> Rk =
            Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>::Identity(dim, dim);

        Rk(k,   k)   =  c;
        Rk(k,   k+1) = -s;
        Rk(k+1, k)   =  s;
        Rk(k+1, k+1) =  c;

        R = R * Rk;   // compose rotations
    }

    // Extract A and b
    auto A = P.get_mat();
    auto b = P.get_vec();

    // Apply A' = A R
    Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> A_rot = A * R;

    return HPOLYTOPE(dim, A_rot, b);
}


//Used to print the ineqialities (if needed)
void print_hpoly(const HPOLYTOPE& P) {
    const auto& A = P.get_mat();   // Matrix: rows = constraints, cols = dimension
    const auto& b = P.get_vec();   // Vector: one entry per inequality

    unsigned int m = A.rows();
    unsigned int d = A.cols();

    std::cout << "H-representation: A x <= b\n";
    std::cout << "Number of inequalities: " << m << "\n";
    std::cout << "Dimension: " << d << "\n\n";

    for (unsigned int i = 0; i < m; ++i) {
        std::cout << i << ":  ";
        for (unsigned int j = 0; j < d; ++j) {
            std::cout << A(i, j) << " * x" << j;
            if (j < d - 1) std::cout << " + ";
        }
        std::cout << " <= " << b(i) << "\n";
    }
    std::cout << std::endl;
}


// Used to write final samples to a file
void write_to_file(std::string filename, std::vector<Point> const& randPoints) {
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); //save and redirect
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf); //reset to standard output again
}

// walk policy
PushBackWalkPolicy push_back_policy;

//Useful class to count time and avoid repeating auto chrono etc
//Add a label when calling an object generator_timer.stop(""); to see time at every loop
class Timer {
public:
    Timer(const std::string& name = "") : walk_name(name), total_time(0.0) {}

    void start() { start_time = std::chrono::steady_clock::now(); }

    // Stop and accumulate
    double stop(const std::string& label = "") {
        auto end_time = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(end_time - start_time).count();
        total_time += elapsed;
        if (!label.empty())
            std::cout << "[" << walk_name << "] " << label << " = " << elapsed << " s\n";
        return elapsed;
    }

    double get_total_time() const { return total_time; }

private:
    std::string walk_name;
    std::chrono::steady_clock::time_point start_time;
    double total_time;
};

// Function to dynamically compute the batch size based on the walk type and dimension
unsigned int compute_batch_size(const std::string& walk_name, unsigned int dim) {

    // Default base size
    unsigned int batch_size = 1000;

    // Adaptive logic by walk type
    if (walk_name == "BallWalk") {
        //Cube
        //batch_size = dim * 150 + 1500;
        //Simplex
        //batch_size = dim * 500 + 1500;
        //Birkhoff
        batch_size = dim * 1000 + 1500;
    }
    else if (walk_name == "AcceleratedBilliardWalk") {
        //Cube
        //batch_size = dim * 2 + 200;
        //Simplex
        //batch_size = dim * 1 + 1500;
        //Birkhoff
        batch_size = dim * 4 + 2000;
    }
    else if (walk_name == "BilliardWalk") {
        //Cube
        //batch_size = dim * 4 + 300;
        //Simplex
        //batch_size = dim * 2 + 3000;
        //Birkhoff
        batch_size = dim * 20 + 2500;
    }
    else if (walk_name == "SparseBilliardWalk") {
        //Cube
        //batch_size = dim * 2 + 200;
        //Simplex
        //batch_size = dim * 1 + 1500;
        //Birkhoff
        batch_size = dim * 2 + 2000;
    }
    else if (walk_name == "CDHRWalk") {
        //Cube
        //batch_size = dim * 1 + 500;
        //Simplex
        //batch_size = dim * 2 + 2400;
        //Birkhoff
        batch_size = dim * 3 + 2400;
        //Biology
        //batch_size = dim * 100 + 2400;
     }
    else if (walk_name == "RDHRWalk") {
        //Cube
        //batch_size = dim * 15 + 500;
        //Simplex
        //batch_size = dim * 350 + 1500;
        //Birkhoff
        batch_size = dim * 350 + 1500;
     } 
    else if (walk_name == "DikinWalk") {
        batch_size = dim * 100 + 1000;
    }
    else if (walk_name == "JohnWalk") {
        batch_size = dim * 250 + 1000;
    }
    else if (walk_name == "VaidyaWalk") {
        batch_size = dim * 150 + 1000;
    }     
       

    return batch_size;
}

// Converts a vector of Points into an Eigen matrix
template <typename MT>
MT vector_to_eigen(const std::vector<Point>& someSamples) {
    MT samples(someSamples[0].dimension(), someSamples.size());
    for (unsigned int jj = 0; jj < someSamples.size(); ++jj)
        samples.col(jj) = someSamples[jj].getCoefficients();
    return samples;
}

// Computes ESS on a given Eigen matrix and returns the min ESS
template <typename NT, typename VT, typename MT>
unsigned int compute_ess(const MT& samples) {
    unsigned int min_ess = 0;
    VT ess_vector = effective_sample_size<NT, VT, MT>(samples, min_ess);
    //std::cout << "Current ESS is: " << min_ess << "\n";
    return min_ess;
}

// Computes PSRF for all accumulated points (no timing inside)
template <typename NT, typename VT, typename MT>
double compute_psrf(const std::vector<Point>& someSamples) {
    // Convert to Eigen matrix
    MT finalSamples(someSamples[0].dimension(), someSamples.size());
    for (unsigned int jj = 0; jj < someSamples.size(); ++jj)
        finalSamples.col(jj) = someSamples[jj].getCoefficients();

    // Compute PSRF
    VT psrf = univariate_psrf<NT, VT, MT>(finalSamples);
    double max_psrf = psrf.maxCoeff();

    return max_psrf;
}

// We choose the walk_len appropriate for each method
unsigned int set_walk_len(const std::string& walk_name, unsigned int dim) {
    static const std::unordered_map<std::string, std::function<unsigned int(unsigned int)>> walk_len_map = {
        {"BallWalk", [](unsigned int dim){ return dim*2; }},
        {"BilliardWalk", [](unsigned int){ return 1; }},
        {"AcceleratedBilliardWalk", [](unsigned int){ return 1; }},
        {"CDHRWalk", [](unsigned int dim){ return dim*2; }},
        {"RDHRWalk", [](unsigned int dim){ return dim*2; }},
        {"DikinWalk", [](unsigned int dim){ return 50+dim; }},
        {"JohnWalk", [](unsigned int dim){ return 50+dim; }},
        {"VaidyaWalk", [](unsigned int dim){ return 50+dim; }}
    };

    if (auto it = walk_len_map.find(walk_name); it != walk_len_map.end())
        return it->second(dim);

    return 1; // default fallback
}

template <typename WalkType>
bool sample_using_walk(HPOLYTOPE& Polytope, 
                       Point const& start_point, 
                       RNGType& rng, 
                       unsigned int target_ESS,
                       const std::string& walk_name,
                       double time_limit_sec) 
{
    Point starting_point = start_point; 
    unsigned int dim = Polytope.dimension();

    Timer t;
    //Uncommenct these 2 lines and comment lines 269-270 and 334-349 to swap back to manual batch size
    unsigned int batch_size = compute_batch_size(walk_name, dim);
    std::cout << "\n[" << walk_name << "] Current batch_size: " << batch_size << "\n";
    //unsigned int batch_size = target_ESS*5;
    //std::cout << "\n[" << walk_name << "] Using dynamic batch size. Initial batch_size: " << batch_size << "\n";

    unsigned int walk_len = set_walk_len(walk_name, dim);
    std::cout << "[" << walk_name << "] Current walk_len: " << walk_len << "\n";
    
    unsigned int current_ESS = 0;
    unsigned int loop_step = 1;
    bool failed_to_converge = false;
    bool timed_out = false; 

    Timer eigen_timer(walk_name);
    Timer ess_timer(walk_name);
    Timer generator_timer(walk_name);
    Timer psrf_timer(walk_name);

    typedef RandomPointGenerator<WalkType> Generator;
    std::vector<Point> allSamples;
    allSamples.reserve(target_ESS * 10); 

    while (current_ESS < target_ESS) {
        
        // CHECK TIME LIMIT
        if (generator_timer.get_total_time() > time_limit_sec) {
            std::cout << "[" << walk_name << "] TIMEOUT (" 
                      << generator_timer.get_total_time() << "s > " 
                      << time_limit_sec << "s). Stopping."
                      << std::string(20, ' ') << "\n"; // Adds 20 spaces to clear the line;
            timed_out = true;
            failed_to_converge = true;
            break;
        }

        //unsigned int num_points = batch_size;
        std::vector<Point> batchPoints;
        
        generator_timer.start();
        Generator::apply(Polytope, starting_point, batch_size, walk_len, batchPoints, push_back_policy, rng);
        generator_timer.stop(""); // Accumulates time

        if (!batchPoints.empty()) {
            starting_point = batchPoints.back();
        }
        allSamples.insert(allSamples.end(), batchPoints.begin(), batchPoints.end());

        // Check ESS
        if (allSamples.size() >= target_ESS) {
            eigen_timer.start();
            MT samples = vector_to_eigen<MT>(allSamples);
            eigen_timer.stop("");
            
            ess_timer.start();
            current_ESS = compute_ess<NT, VT, MT>(samples);
            ess_timer.stop("");

            if (current_ESS >= target_ESS) break;
            
            // Optional: Print progress
            std::cout << "[" << walk_name << "] Samples: " << allSamples.size() 
                      << " | ESS: " << current_ESS 
                      << " | Time: " << generator_timer.get_total_time() << "s\r" << std::flush;
        } else {
            current_ESS = 0;
        }

        // // // THE DYNAMIC batch size logic
        // // Calculate how useful each sample is
        // double ess_per_sample = (double)current_ESS / (double)allSamples.size();
        
        // if (ess_per_sample > 1e-6) {
        //     unsigned int remaining_ESS = target_ESS - current_ESS;
        //     // Predict needed samples + 10% safety buffer
        //     unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 1.1);
            
        //     // Don't let the batch size get smaller than the initial guess, 
        //     // but allow it to scale up to finish the job in the next loop.
        //     batch_size = std::max(batch_size, samples_needed);
        // } else {
        //     // If efficiency is essentially zero, double the batch to try and find signal
        //     batch_size *= 2; 
        // }

        // Failsafe for infinite loops
        if (loop_step > 50 && current_ESS < target_ESS) {
            std::cout << "\n[" << walk_name << "] I am sorry. I did not converge!\n";
            failed_to_converge = true;
            break;
        }
        loop_step++;
    }
    std::cout << "[" << walk_name << "] Samples: " << allSamples.size() 
              << " | ESS: " << current_ESS 
              << " | Time: " << generator_timer.get_total_time() << "s\r" << std::flush;

    // PRINT RESULTS TO CONSOLE
    if(!failed_to_converge) {
        std::cout << "\n[" << walk_name << "] DONE. Final ESS: " << current_ESS << "\n";
        std::cout << "[" << walk_name << "] Total generation time: " << generator_timer.get_total_time() << " s\n";
    }

    // SAVE TO FILE 
    // Format: Dimension, WalkName, Time, Points, ESS
    std::ofstream outfile;
    outfile.open("benchmark_results.txt", std::ios_base::app); // Append mode
    
    if (outfile.is_open()) {
        outfile << Polytope.dimension() << ", " 
                << walk_name << ", "
                << generator_timer.get_total_time() << ", "
                << allSamples.size() << ", "
                << current_ESS << "\n";
        outfile.close();
        std::cout << "[" << walk_name << "] Results saved to benchmark_results.txt"
                  << std::string(20, ' ') << "\n"; // Adds 20 spaces to clear the line;
    } else {
        std::cerr << "Unable to open file to write results!\n";
    }

    // PSRF calculation
    psrf_timer.start();
    double max_psrf = compute_psrf<NT, VT, MT>(allSamples);
    psrf_timer.stop("Total PSRF time");

    // PSRF result
    std::cout << "[" << walk_name << "] Total time to calculate ESS = " << ess_timer.get_total_time() << " s\n";
    std::cout << "[" << walk_name << "] PSRF = " << max_psrf << std::endl;

    // K-S statistical test ////////////////////////////////////////////////////
    MT samples_mat = vector_to_eigen<MT>(allSamples);
    int computed_thin = static_cast<int>(samples_mat.cols() / current_ESS);
    int thin_factor = std::max(10, computed_thin * 2);
    if (thin_factor < 1) thin_factor = 1;

    auto [ks_stat, p_val, observed, expected] = global_scaling_test(Polytope, samples_mat, thin_factor);  

    std::cout << "[" << walk_name << "] KS Statistic: " << ks_stat << "\n";
    std::cout << "[" << walk_name << "] P-Value:      " << p_val << "\n"; 

    // Comment-uncomment to see full volume analysis
    // Print the Shell Analysis Table (using the unpacked variables)
    // std::cout << "Volume Shells Analysis (Expected vs Observed):\n";
    // std::cout << "Exp Vol% | Obs Vol% | Deviation\n";
    // std::cout << "-------------------------------\n";

    // double max_dev_perc = 0.0;
    // for(size_t i = 0; i < 10; ++i) {
    //    double dev = (observed[i] - expected[i]) * 100.0;
        
    //    if(std::abs(dev) > std::abs(max_dev_perc)) max_dev_perc = dev;
        
    //    printf("  %4.1f%%  |  %4.1f%%  | %+6.2f%%\n", 
    //           expected[i]*100.0, observed[i]*100.0, dev);
    // }
    // std::cout << "-------------------------------\n";

    if (timed_out || failed_to_converge) return false; 

    // Uncomment to write samples to file
    //write_to_file("All_Samples.txt", allSamples);
    //std::cout << "All samples written to ALL_Samples.txt file. DONE" << std::endl;

    return true;
}

///////////////RIEMMANIAN//////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename PolytopeType, typename RNGType>
bool sample_using_crhmc(PolytopeType& HP, 
                        typename PolytopeType::PointType& /*center*/, 
                        RNGType& rng, 
                        unsigned int target_ESS,
                        const std::string& walk_name) 
{
    using NT = double;
    using Point = typename PolytopeType::PointType;
    using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
    using MT = typename PolytopeType::MT;
    
    using Func = ZeroScalarFunctor<Point>;
    using Grad = ZeroFunctor<Point>;
    using Hess = ZeroFunctor<Point>;

    // --- Configuration ---
    double current_efficiency = 1.0 / 10.0; 
    int n_burns = 1000; 
    int walk_len = 1; // Thinning
    
    // Global container for ALL samples across batches
    std::vector<Point> all_samples;
    all_samples.reserve(target_ESS * 20);

    Timer total_timer(walk_name);
    total_timer.start();

    double current_ESS = 0.0;
    int batch_count = 1;

    std::cout << "[" << walk_name << "] Target ESS: " << target_ESS << "\n";

    // --- The Smart Loop ---
    while (current_ESS < target_ESS) {
        
        // 1. Calculate how many samples we need
        double missing_ESS = target_ESS - current_ESS;
        
        // "Smart Batching": Estimate samples needed based on current efficiency
        // We add a 10% buffer (1.1) to try and finish in this batch
        int n_samples_needed = static_cast<int>((missing_ESS / current_efficiency) * 1.1);
        
        // Safety clamps
        if (n_samples_needed < 1000) n_samples_needed = 1000; 
        if (n_samples_needed > 100000) n_samples_needed = 100000; // Cap to prevent memory explosion

        std::cout << "[" << walk_name << "][Batch " << batch_count << "] Requesting " 
                  << n_samples_needed << " samples (Efficiency: " << current_efficiency << ")\n";

        // Setup Helper Objects (Must be fresh per run)
        Func* f = new Func;
        Grad* g = new Grad;
        std::list<Point> batch_list;

        // EXECUTE CRHMC
        // Note: passing 'rng' ensures the random sequence continues, making this valid
        execute_crhmc<PolytopeType, RNGType, std::list<Point>, Grad, Func, Hess, CRHMCWalk, 1>(
            HP, rng, batch_list, walk_len, n_samples_needed, n_burns, g, f
        );

        delete f;
        delete g;

        // Merge Samples
        // We move elements from list to vector to avoid copying
        all_samples.insert(all_samples.end(), batch_list.begin(), batch_list.end());

        // 5. Update Statistics
        // We must convert ALL samples to matrix to calculate total ESS
        MT samples_matrix = MT(HP.dimension(), all_samples.size());
        for (size_t i = 0; i < all_samples.size(); ++i) {
            samples_matrix.col(i) = all_samples[i].getCoefficients();
        }

        // Calculate ESS
        current_ESS = compute_ess<NT, VT, MT>(samples_matrix);

        // Update Efficiency for next loop
        // Efficiency = ESS / Total_Raw_Samples
        if (all_samples.size() > 0) {
            current_efficiency = current_ESS / all_samples.size();
        }

        std::cout << "[" << walk_name << "][Batch " << batch_count << "] Current ESS: " << current_ESS 
                  << " / " << target_ESS << "\n";

        // Break if we are stuck (Efficiency drops too low)
        if (current_efficiency < 0.0001 && all_samples.size() > 5000) {
            std::cout << "[" << walk_name << "] CRITICAL: Efficiency too low. Stopping.\n";
            break;
        }

        batch_count++;
    }
    
    total_timer.stop("");

    // --- Final Reporting & Tests ---
    
    // 1. PSRF
    Timer psrf_timer(walk_name);
    psrf_timer.start();
    double max_psrf = compute_psrf<NT, VT, MT>(all_samples);
    psrf_timer.stop("");

    // 2. KS Test
    MT final_matrix = MT(HP.dimension(), all_samples.size());
    for (size_t i = 0; i < all_samples.size(); ++i) {
        final_matrix.col(i) = all_samples[i].getCoefficients();
    }
    
    // KS Logic
    double safe_ess = (current_ESS > 0.0) ? current_ESS : 1.0;
    int computed_thin = static_cast<int>(final_matrix.cols() / safe_ess);
    int thin_factor = std::max(10, computed_thin * 2);
    if (thin_factor < 1) thin_factor = 1;

    auto [ks_stat, p_val, observed, expected] = global_scaling_test(HP, final_matrix, thin_factor);

    std::cout << "------------------------------------------------\n";
    std::cout << "[" << walk_name << "] Total time     : " << total_timer.get_total_time() << " s\n";
    std::cout << "[" << walk_name << "] Total Batches  : " << (batch_count - 1) << "\n";
    std::cout << "[" << walk_name << "] Total Samples  : " << all_samples.size() << "\n";
    std::cout << "[" << walk_name << "] Final ESS      : " << current_ESS << "\n";
    std::cout << "[" << walk_name << "] Final PSRF     : " << max_psrf << "\n";
    std::cout << "[" << walk_name << "] KS Statistic   : " << ks_stat << "\n";
    std::cout << "[" << walk_name << "] P-Value        : " << p_val << "\n";
    std::cout << "------------------------------------------------\n";
    
    // Uncomment to write samples to file
    //write_to_file("All_Samples.txt", all_samples);
    //std::cout << "All samples written to ALL_Samples.txt file. DONE" << std::endl;
    return true;
}
/////////////////End of Riemannian///////////////////////////////////////////////////////////////////////////////////////////////////////

typedef BallWalk::template Walk<HPOLYTOPE, RNGType> BallWalkType;
typedef BilliardWalk::template Walk<HPOLYTOPE, RNGType> BilliardWalkType;
typedef AcceleratedBilliardWalk::template Walk<HPOLYTOPE, RNGType> AcceleratedBilliardWalkType;
typedef CDHRWalk::template Walk<HPOLYTOPE, RNGType> CDHRWalkType;
typedef DikinWalk::template Walk<HPOLYTOPE, RNGType> DikinWalkType;
typedef JohnWalk::template Walk<HPOLYTOPE, RNGType> JohnWalkType;
typedef RDHRWalk::template Walk<HPOLYTOPE, RNGType> RDHRWalkType;
typedef VaidyaWalk::template Walk<HPOLYTOPE, RNGType> VaidyaWalkType;
typedef SparseBilliardWalk::template Walk<HPOLYTOPE, RNGType> SparseBilliardaWalkType;

int main(int argc, char const *argv[]) {

    // You can adjust dimensions as needed
    std::vector<unsigned int> dimensions = {500};
    // Select target ESS
    unsigned int target_ESS = 800;

    //seed for order polytopes
    int base_seed = 42;
    
    // Time Limit: 1 Hour (3600 seconds)
    // If a method exceeds this, it stops and is skipped for all larger dimensions.
    double TIME_LIMIT_SEC = 1200.0; 

    // Choose what methods to use. True is used, false is not used.
    std::map<std::string, bool> active_methods;
    active_methods["AcceleratedBilliardWalk"] = true;
    active_methods["CDHRWalk"]                = true;
    active_methods["BallWalk"]                = false;
    active_methods["BilliardWalk"]            = false;
    active_methods["RDHRWalk"]                = false;
    active_methods["CRHMCWalk"]               = false;
    active_methods["SparseBilliardWalk"]      = false; 

    // Initialize output file (overwrite old one)
    std::ofstream outfile("benchmark_results.txt");
    outfile << "Dim, Method, Time(s), Points, ESS\n";
    outfile.close();

    std::cout << "Starting Benchmark.\n"; 
    std::cout << "Target ESS: " << target_ESS << "\n";
    std::cout << "Time Limit: " << TIME_LIMIT_SEC << "s per method.\n";

    //////CUSTOM Polytopes in A*x<=b form/////////////*******************************************////////////////////////////////////////
    // TOGGLE THIS: Set to 'true' for your custom CSVs, 'false' for the Cube/Simplex/... benchmark
    bool USE_CUSTOM_MODEL = false; 

    HPOLYTOPE custom_polytope; // Placeholder for the loaded model

    if (USE_CUSTOM_MODEL) {
        try {
            // Load the model ONCE before the loop. Place the csv files in the build folder.
            custom_polytope = load_custom_polytope<HPOLYTOPE>("agg_A.csv", "agg_b.csv");
            
            // Overwrite dimensions list to run exactly ONCE for the model's dimension
            dimensions = { static_cast<unsigned int>(custom_polytope.dimension()) };
            
            std::cout << ">>> CUSTOM MODE ACTIVATED: Loaded model with " << dimensions[0] << " dimensions.\n";
        } catch (const std::exception& e) {
            std::cerr << "CRITICAL ERROR: Could not load custom files. " << e.what() << std::endl;
            return 1;
        }
    } 
    //////End of Custom Polytopes///////////////////*******************************************////////////////////////////////////////

    for (auto dim : dimensions) {
        
        std::cout << "\n" << std::string(40, '=') << "\n";
        std::cout << "*** Running for dimension " << dim << " ***\n";

        HPOLYTOPE Polytope_simple;
        HPOLYTOPE Polytope;
        if (USE_CUSTOM_MODEL) {
            Polytope = custom_polytope;
        } else {
            // Generate desired Polytope 
            //Polytope_simple = generate_cube<HPOLYTOPE>(dim, false);
            //Polytope_simple = generate_birkhoff<HPOLYTOPE>(dim);
            //Polytope_simple = generate_cross<HPOLYTOPE>(dim, false);
            //Polytope_simple = generate_skinny_cube<HPOLYTOPE>(dim,false);
            //Polytope_simple = generate_simplex<HPOLYTOPE>(dim, false);

            // Generate order polytopes
            unsigned int m = 3 * dim;
            int current_seed = base_seed + dim;
            std::cout << "\nCreating order polytope...\n";
            Polytope_simple = random_orderpoly<HPOLYTOPE, double>(dim, m, current_seed);
            
            double angle = 53.0 * M_PI / 180.0; //53 deg
            //double angle = 0.0 * M_PI / 180.0;  //0 deg
            Polytope = rotate_all_dims(Polytope_simple, angle);

            //print_hpoly(Polytope); //Uncomment to print the polytope
        }
        // Setup RNG and Starting Point
        RNGType rng(Polytope.dimension());
        auto inner = Polytope.ComputeInnerBall();
        Point center = inner.first;

        //Use the following lines to force the center (starting point) be 0
        // point<Cartesian<double>> center(dim);
        // center.set_to_origin();
        // std::cout << "Starting point: ";
        // center.print();

        // =========================================================
        // ACCELERATED BILLIARD WALK
        // =========================================================
        if (active_methods["AcceleratedBilliardWalk"]) {
            bool success = sample_using_walk<AcceleratedBilliardWalkType>(
                Polytope, center, rng, target_ESS, "AcceleratedBilliardWalk", TIME_LIMIT_SEC
            );
            if (!success) {
                active_methods["AcceleratedBilliardWalk"] = false;
                std::cout << "!!! Disabling AcceleratedBilliardWalk for future dimensions.\n";
            }
        } else {
            std::cout << "\n[AcceleratedBilliardWalk] Skipping (previously timed out).\n";
        }

        // =========================================================
        // CDHR (Coordinate Directions Hit-and-Run)
        // =========================================================
        if (active_methods["CDHRWalk"]) {
            bool success = sample_using_walk<CDHRWalkType>(
                Polytope, center, rng, target_ESS, "CDHRWalk", TIME_LIMIT_SEC
            );
            if (!success) {
                active_methods["CDHRWalk"] = false;
                std::cout << "!!! Disabling CDHR for future dimensions.\n";
            }
        } else {
            std::cout << "\n[CDHRWalk] Skipping (previously timed out).\n";
        }

        // =========================================================
        // BILLIARD WALK
        // =========================================================
        if (active_methods["BilliardWalk"]) {
            bool success = sample_using_walk<BilliardWalkType>(
                Polytope, center, rng, target_ESS, "BilliardWalk", TIME_LIMIT_SEC
            );
            if (!success) {
                active_methods["BilliardWalk"] = false;
                std::cout << "!!! Disabling Billiard Walk for future dimensions.\n";
            }
        } else {
            std::cout << "\n[BilliardWalk] Skipping (previously timed out).\n";
        }

        // =========================================================
        // RDHR (Random Directions Hit-and-Run)
        // =========================================================
        if (active_methods["RDHRWalk"]) {
            bool success = sample_using_walk<RDHRWalkType>(
                Polytope, center, rng, target_ESS, "RDHRWalk", TIME_LIMIT_SEC
            );
            if (!success) {
                active_methods["RDHRWalk"] = false;
                std::cout << "!!! Disabling RDHR for future dimensions.\n";
            }
        } else {
            std::cout << "\n[RDHRWalk] Skipping (previously timed out).\n";
        }

        // =========================================================
        // BALL WALK
        // =========================================================
        if (active_methods["BallWalk"]) {
            bool success = sample_using_walk<BallWalkType>(
                Polytope, center, rng, target_ESS, "BallWalk", TIME_LIMIT_SEC
            );
            if (!success) {
                active_methods["BallWalk"] = false;
                std::cout << "!!! Disabling Ball Walk for future dimensions.\n";
            }
        } else {
            std::cout << "\n[BallWalk] Skipping (previously timed out).\n";
        }

        // =========================================================
        // SparseBilliardWalk
        // =========================================================
        if (active_methods["SparseBilliardWalk"]) {
            bool success = sample_using_walk<SparseBilliardaWalkType>(
                Polytope, center, rng, target_ESS, "SparseBilliardWalk", TIME_LIMIT_SEC
            );
            if (!success) {
                active_methods["SparseBilliardWalk"] = false;
                std::cout << "!!! Disabling Sparse Billiard Walk for future dimensions.\n";
            }
        } else {
            std::cout << "\n[SparseBilliardWalk] Skipping (previously timed out).\n";
        }

        // =========================================================
        // CRHMC WALK
        // =========================================================
        if (active_methods["CRHMCWalk"]) {
            std::cout << "Starting CRHMC..." << std::endl;
            bool success = sample_using_crhmc(Polytope, center, rng, target_ESS, "CRHMCWalk");
            if (!success) {
                active_methods["CRHMCWalk"] = false;
                std::cout << "!!! Disabling CRHMCWalk for future dimensions.\n";
            }
        } else {
            std::cout << "\n[CRHMCWalk] Skipping (previously timed out).\n";
        }

    } // End dimension loop

    std::cout << "\nBenchmark Complete.\n";
    return 0;
}