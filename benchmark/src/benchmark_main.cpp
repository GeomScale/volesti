#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

#include "../include/core_types.hpp"
#include "../include/walk_parameters.hpp"
#include "../include/benchmark_utils.hpp"
#include "../include/geometry_utils.hpp"
#include "../include/walk_result.hpp"
#include "../include/walk_registry.hpp"
#include "../include/menu.hpp"

#include "../include/polytope_generation.hpp"
#include "known_polytope_generators.h"
#include "order_polytope_generator.h"
#include "inscribed_ellipsoid_rounding.hpp"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char** argv) {

    initialize_all_walks();

    string config_file = "../config/walk_config.json"; // Base default
    unsigned int dimension;
    string walk_choice;
    string polytope_cli_choice;

    // We first pre-parse only the config file flag so we can load the JSON defaults first
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            config_file = argv[i + 1];
            break;
        }
    }

    // Load the Configuration
    cout << "Loading configuration from: " << config_file << "\n";
    BenchmarkConfig config = load_benchmark_config(config_file);

    // Setup Command Line Options (using JSON values as our defaults)
    po::options_description desc("Benchmark Options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("config,c", po::value<string>(&config_file)->default_value(config_file), "Path to JSON config")
        ("dim,d", po::value<unsigned int>(&dimension)->default_value(config.dimension), "Dimension of the polytope")
        ("polytope,p", po::value<string>(&polytope_cli_choice)->default_value(config.polytope_choice), "Choose polytope: Cube, Simplex, Birkhoff, Cross, OrderPolytope, Custom")
        ("walk,w", po::value<string>(&walk_choice)->default_value("All"), 
            "Specific walk to run, or 'All'. Valid options:\n"
            "  - BallWalk\n"
            "  - BilliardWalk\n"
            "  - AcceleratedBilliardWalk\n"
            "  - SparseBilliardWalk\n"
            "  - CDHRWalk\n"
            "  - RDHRWalk\n"
            "  - DikinWalk\n"
            "  - JohnWalk\n"
            "  - VaidyaWalk\n"
            "  - GaussianBallWalk\n"
            "  - GaussianCDHRWalk\n"
            "  - BilliardShakeAndBakeWalk\n"
            "  - ShakeAndBakeWalk\n"
            "  - BCDHRWalk\n"
            "  - BRDHRWalk\n"
            "  - CRHMCWalk");

    // Error-Checking Block
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const exception& e) {
        cerr << "Error parsing arguments: " << e.what() << "\n";
        return 1;
    }

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }

    config.polytope_choice = polytope_cli_choice;

    // --- NEW MENU LOGIC ---
    if (config.show_menu) {
        bool continue_to_benchmark = run_interactive_menu(config, walk_choice);
        if (!continue_to_benchmark) {
            return 0; // Exit gracefully if they chose 3
        }
    }
    // ---------------------------

    std::vector<unsigned int> dimensions_to_run;
    if (config.polytope_choice == "Custom") {
        dimensions_to_run = { 0 };
    } else {
        dimensions_to_run = config.dimensions;
    }

    for (unsigned int current_dim : dimensions_to_run) {

        config.dimension = current_dim;
        HPOLYTOPE Polytope_simple;
        try {
            Polytope_simple = create_polytope(config.polytope_choice, config.dimension, config);
        } 
        catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            std::cerr << ">>> Error. Please fix the config or file paths.\n";
            return 1; 
        }

        config.dimension = Polytope_simple.dimension();
        dimension = config.dimension; 

        // Print basic info
        cout << "Target ESS: " << config.target_ESS << "\n";
        cout << "Dimension: " << dimension << "\n";
        cout << "Polytope: " << config.polytope_choice << "\n";

        double angle = config.angle; 
        cout << "Rotation angle is: " << angle << "\n";
        cout << "Dynamic batch size is on: " << config.use_dynamic_batch << "\n";
        cout << "Rounding is on: " << config.rounding << "\n";
        cout << "Auto-walk is on: " << (config.auto_walk ? "true" : "false") << "\n";

        cout << "\n" << string(40, '=') << "\n";
        cout << "*** Running for dimension " << dimension << " ***\n";
        
        // We need this copy to pass the original polytope to the metrics if roundeing was on.
        HPOLYTOPE Polytope_rotated = rotate_all_dims(Polytope_simple, angle);
        HPOLYTOPE Polytope = Polytope_rotated;
        
        auto inner = Polytope.ComputeInnerBall();
        Point center = inner.first;

        // Variables to store transformation data
        MT T;
        VT shift;
        NT round_val = 1.0;

        // ****ROUNDING***** 
        if (config.rounding) {
            cout << "[ROUNDING] Rounding is enabled. Applying max inscribed ellipsoid rounding...\n";
            
            // Pass the pre-computed center into the rounding function
            // Use john ellispoid
            // Polytope is passed by reference so we shouldnt need to define a new one.
            auto rounding_result = inscribed_ellipsoid_rounding<MT, VT, NT>(Polytope, center);
            
            // Unpack the transformation data
            T = std::get<0>(rounding_result);
            shift = std::get<1>(rounding_result);
            round_val = std::get<2>(rounding_result);
            
            // Since rounding shifts the polytope to the origin we will use 0,0,0,0,0 ... as center.
            center = Point(VT::Zero(Polytope.dimension()));
            
            cout << "[ROUNDING] Rounding complete. Round value: " << round_val << "\n\n";
        }
        // --------------------------

        // Setup RNG 
        RNGType rng(Polytope.dimension());

        auto run_method = [&](const string& method_name) {

            // Access registry
            auto& registry = get_walk_registry();
            auto walk_it = registry.find(method_name);

            if (walk_it != registry.end()) {
                
                WalkResult result = walk_it->second(
                    Polytope,
                    center,
                    rng,
                    config,
                    method_name
                );

                if (!result.samples.empty()) {

                    // Reverse rounding
                    if (config.rounding) {
                        for (auto& pt : result.samples) {
                    
                            // Apply the reverse transformation: T * vector + shift
                            pt = T * pt.getCoefficients() + shift; 
                        }
                    }
                    // ----------------------------------

                    // Write samples to txt file for later use
                    if (config.write_to_file) {
                        std::filesystem::create_directory("results");
                        std::string filename = "results/" + config.polytope_choice + "_" + 
                                               std::to_string(dimension) + "_" + 
                                               method_name + "_samples.txt";
                        
                        std::cout << "[" << method_name << "] Saving " << result.samples.size() 
                                << " points to " << filename << "...\n";
                                
                        write_to_file(filename, result.samples);
                        
                        std::cout << "[" << method_name << "] File saved successfully.\n";
                    }

                    // Process results
                    process_and_print_results(
                        result.samples, 
                        Polytope_rotated, 
                        method_name, 
                        result.generation_time, 
                        result.final_ess,
                        result.ess_time,
                        result.walk_len        
                    );
                } else {
                    cout << "!!! " << method_name << " failed to generate points.\n";
                }
                
            } else {
                cout << "!!! Unknown walk type skipped: " << method_name << "\n";
                return;
            }
        };
        // Auto-walk handles the selection if enabled
        if (config.auto_walk) {
            string auto_selected_walk = determine_auto_walk(dimension);
            cout << "\n[Auto-Walk] Dimension " << dimension << " overriding config to run: " << auto_selected_walk << "\n";
            run_method(auto_selected_walk);
        }
        // If the user picked "All", iterate through the JSON keys. Otherwise, just run the one they requested.
        else if (walk_choice == "All") {
            for (const auto& pair : config.walk_settings) {

                if(pair.second.enabled) {
                    run_method(pair.first);
                } 
                else {
                // cout << "--- Skipping " << pair.first << " (Disabled in JSON) ---\n";
                }
            }
        } else {
            run_method(walk_choice);
        }
    }
    cout << "\nBenchmark Complete.\n";
    return 0;
}