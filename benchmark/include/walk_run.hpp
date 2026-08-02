#pragma once

#include "core_types.hpp"
#include "walk_parameters.hpp"
#include "benchmark_utils.hpp" 
#include "walk_result.hpp"          
#include "walk_adapters.hpp"
#include "diagnostics.hpp"
#include "progress_bar.hpp"
#include "dynamic_batch_size.hpp"

#include "sampling/random_point_generators.hpp"
#include <iostream>
#include <filesystem>

// In this file you can find sample using walk function used to generate
// points from a polytope using a specified random walk algorithm.

template <typename WalkType>
WalkResult sample_using_walk(HPOLYTOPE& Polytope,  
                             Point const& start_point, 
                             RNGType& rng, 
                             const BenchmarkConfig& config,
                             const std::string& walk_name)
{
    // Initial point and dimension
    Point starting_point = start_point; 
    unsigned int dim = Polytope.dimension();

    // ESS placeholder and loop number
    unsigned int current_ESS = 0;
    unsigned int loop_step = 1;
    unsigned int previous_ESS = 0;
    double estimated_remaining_seconds = -1.0;

    // Random generator and placeholder matrix for samples
    typedef RandomPointGenerator<WalkType> Generator;
    std::vector<Point> allSamples;
    allSamples.reserve(config.target_ESS * 10); 

    // Find initial batch size and walk_len based on config file and user choices
    unsigned int batch_size = get_initial_batch_size(walk_name, config);
    unsigned int walk_len = compute_dynamic_walk_len(walk_name, dim, config);

    // initiallize the walk
    auto walk = WalkAdapter<WalkType>::init(Polytope, starting_point, config, rng);

    // Timer objects that stored the time
    Timer walk_timer(walk_name);
    Timer ess_timer(walk_name);
    walk_timer.start();
    Timer global_timer("TotalTime");
    global_timer.start();

    // Main while loop. We sample intil we hit target ESS.
    while (current_ESS < config.target_ESS && global_timer.get_total_time() < config.time_limit_sec) {
        
        std::vector<Point> batchPoints;

        // here we calculate the live time estimate for the upcoming batch 
        double elapsed_walk_time = walk_timer.get_total_time();
        unsigned int samples_so_far = allSamples.size();

        if (samples_so_far > 0 && elapsed_walk_time > 0.01 && current_ESS > 0) {
            // how many samples we generate per second
            double samples_per_sec = static_cast<double>(samples_so_far) / elapsed_walk_time;
            // how many samples are still needed to reach Target ESS
            double remaining_ess_needed = static_cast<double>(config.target_ESS) - current_ESS;
            // Samples needed = (Remaining ESS) * (Average samples required per 1 ESS)
            double estimated_remaining_samples = remaining_ess_needed * (static_cast<double>(samples_so_far) / current_ESS);
            // ETA = Remaining Samples / Sample Rate
            estimated_remaining_seconds = estimated_remaining_samples / samples_per_sec;
        } else {
            estimated_remaining_seconds = -1.0; // Still on the first batch / no ESS data yet
        }
            
        // Progress Bar
        unsigned int chunk_size = 250; // How many samples to generate before updating the bar
        unsigned int generated_this_batch = 0;

        if constexpr (WalkAdapter<WalkType>::supports_chunking) {
                
                unsigned int chunk_size = 250; 
                unsigned int generated_this_batch = 0;
                unsigned int samples_before_this_batch = allSamples.size();

                while (generated_this_batch < batch_size) {
                    unsigned int current_chunk = std::min(chunk_size, batch_size - generated_this_batch);
                    std::vector<Point> chunkPoints;
                    
                    // The main call. The same call is valid for all methods
                    WalkAdapter<WalkType>::apply_batch(
                        walk, Polytope, starting_point, current_chunk, walk_len, chunkPoints, config, rng, walk_timer
                    );

                    if (!chunkPoints.empty()) {
                        starting_point = chunkPoints.back();
                    }
                    allSamples.insert(allSamples.end(), chunkPoints.begin(), chunkPoints.end());
                    generated_this_batch += chunkPoints.size();
                    double global_elapsed_seconds = global_timer.get_total_time();

                    // Draw progress bar with mixing ratio
                    draw_progress_bar(
                        walk_name, 
                        generated_this_batch, 
                        batch_size, 
                        samples_before_this_batch, 
                        current_ESS, 
                        walk_len,
                        global_elapsed_seconds,
                        estimated_remaining_seconds
                    );

                    if (walk_timer.get_total_time() > config.time_limit_sec) {
                        break;
                    }
                }
            } else {
                // RIEMANNIAN METHOD: Massive Single Batch
                std::cout << "\r" << std::string(100, ' ') << "\r[" << walk_name 
                        << "] Generating massive batch of " << batch_size 
                        << " points (Tuning physics engine)..." << std::flush;

                if (estimated_remaining_seconds >= 0.0) {
                    std::cout << " (ETA: " << static_cast<int>(estimated_remaining_seconds) << "s remaining)..." << std::flush;
                }       
                
                std::vector<Point> singleBatchPoints;
                WalkAdapter<WalkType>::apply_batch(
                        walk, Polytope, starting_point, batch_size, walk_len, singleBatchPoints, config, rng, walk_timer
                    );

                if (!singleBatchPoints.empty()) {
                    starting_point = singleBatchPoints.back();
                }
                allSamples.insert(allSamples.end(), singleBatchPoints.begin(), singleBatchPoints.end());
                
                if (walk_timer.get_total_time() > config.time_limit_sec) {
                    std::cout << "\n[" << walk_name << "] Time limit exceeded during batch generation.\n";
                    break;
                }
            }

            std::cout << "\r" << std::string(100, ' ') << "\r[" << walk_name << "] Batch Complete! Calculating ESS..." << std::flush;

            walk_timer.stop("");
            ess_timer.start();

            // Calculate ESS 
            MT samples = vector_to_eigen<MT>(allSamples);
            current_ESS = compute_ess<NT, VT, MT>(samples);

            ess_timer.stop("");
            walk_timer.start();

            std::cout << "\r" << std::string(120, ' ') << "\r[" << walk_name << "] Samples: " << allSamples.size() 
                    << " | ESS: " << current_ESS;

            // If the user has dynamic_batch_size off we end here after "samples" number of samples are generated.
            if (!config.use_dynamic_batch) {
                std::cout << "\n[" << walk_name << "] Generated " << allSamples.size() 
                        << " fixed samples. Stopping without checking ESS.\n";
                break; 
            }

            // Main check. If we passed the target stop immediatelly.
            if (current_ESS >= config.target_ESS) {
                std::cout << "\n[" << walk_name << "] Reached Target ESS. Stopping.\n";
                break;
            }

            if (global_timer.get_total_time() > config.time_limit_sec) {
                std::cout << "\n[" << walk_name << "] Time limit exceeded. Stopping.\n";
                break;
            }

            // THE DYNAMIC BATCH SIZE 
            if (config.use_dynamic_batch) {
                double remaining_time_sec = config.time_limit_sec - global_timer.get_total_time();
                double current_samples_per_sec = (elapsed_walk_time > 0.01) 
                    ? static_cast<double>(allSamples.size()) / elapsed_walk_time 
                    : 0.0;

                batch_size = compute_next_batch_size(
                    config.target_ESS, 
                    current_ESS, 
                    previous_ESS,    
                    allSamples.size(), 
                    batch_size,
                    config.dimension,
                    remaining_time_sec,
                    current_samples_per_sec
                );
            }

            std::cout << " | Next Batch: " << batch_size << std::flush;

            previous_ESS = current_ESS;

            // Failsafe
            if (loop_step > 30 && current_ESS < config.target_ESS) {
                std::cout << "\n[" << walk_name << "] I am sorry. I did not converge!\n";
                break;
            }
            
            loop_step++;
    }

    std::cout << "\n[" << walk_name << "] DONE. Samples Generated: " << allSamples.size() << "\n";
    //global_timer.stop("Total Execution Time");
    double final_gen_time = walk_timer.get_total_time();
    double final_ess_time = ess_timer.get_total_time();
    
    return { allSamples, current_ESS, final_gen_time, final_ess_time, walk_len };
}
