#pragma once

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cstdlib> 
#include <sstream>
#include <algorithm>
#include <cstdlib> 
#include "walk_parameters.hpp"

// Helper function to safely get an integer or quit via 'q'/'Q' from ANY prompt
inline int get_valid_int(const std::string& prompt, int min_val, int max_val) {
    std::string input;
    while (true) {
        std::cout << prompt;
        if (!(std::cin >> input)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "[!] Invalid input. Enter a number or 'q' to quit: ";
            continue;
        }

        // Check for exit command anywhere
        if (input == "q" || input == "Q") {
            std::cout << "\nExiting program. Goodbye!\n";
            std::exit(0);
        }

        // Try parsing the string into an integer safely
        try {
            size_t idx;
            int choice = std::stoi(input, &idx);
            
            // Ensure the entire string was consumed (rejects strings like "5abc")
            if (idx != input.size()) {
                std::cout << "[!] Invalid format. Please enter a valid number or 'q' to quit.\n";
                continue;
            }

            if (choice < min_val || choice > max_val) {
                std::cout << "[!] Out of range. Please choose between " << min_val << " and " << max_val << " (or 'q' to quit).\n";
                continue;
            }
            return choice;
        } catch (const std::exception&) {
            std::cout << "[!] Invalid input. Please enter a valid number or 'q' to quit.\n";
        }
    }
}

// Helper function to get a comma-separated list of integers
inline std::vector<unsigned int> get_valid_dimension_list(const std::string& prompt) {
    std::string input;
    while (true) {
        std::cout << prompt;
        
        // std::ws skips any leftover newlines in the input buffer from previous prompts
        if (!std::getline(std::cin >> std::ws, input)) {
            std::cin.clear();
            std::cout << "[!] Invalid input.\n";
            continue;
        }

        if (input == "q" || input == "Q") {
            std::cout << "\nExiting program. Goodbye!\n";
            std::exit(0);
        }

        std::vector<unsigned int> dims;
        std::stringstream ss(input);
        std::string token;
        bool error = false;

        // Split the string by commas
        while (std::getline(ss, token, ',')) {
            // Strip any accidental spaces (e.g. if they type "10, 20, 30")
            token.erase(std::remove_if(token.begin(), token.end(), ::isspace), token.end());
            
            if (token.empty()) continue;

            try {
                size_t idx;
                int val = std::stoi(token, &idx);
                
                // Ensure it's a valid integer within a reasonable range
                if (idx != token.size() || val < 1 || val > 100000) {
                    error = true;
                    break;
                }
                dims.push_back(static_cast<unsigned int>(val));
            } catch (...) {
                error = true;
                break;
            }
        }

        if (error || dims.empty()) {
            std::cout << "[!] Invalid format. Enter positive integers separated by commas (e.g., 10,20,30).\n";
            continue;
        }

        return dims;
    }
}

inline void show_help() {
    int choice = -1; 

    while (true) {
        // Clear the console screen (cross-platform compatible)
#if defined(_WIN32)
        int ret = std::system("cls");
#else
        int ret = std::system("clear");
#endif
        (void)ret; // Silence the "ignoring return value" warning

        // Print the fixed menu at the top
        std::cout << "\n" << std::string(60, '-') << "\n";
        std::cout << "*************************** HELP ***************************\n";
        std::cout << std::string(60, '-') << "\n";
        std::cout << "Select a topic to learn more about:\n";
        std::cout << "  1. What is an MCMC method?\n";
        std::cout << "  2. What is ESS (Effective Sample Size)?\n";
        std::cout << "  3. What is PSRF?\n";
        std::cout << "  4. What is the KS test?\n";
        std::cout << "  5. Polytope Rounding\n";
        std::cout << "  6. Geometric Methods (Ball, Billiards, Hit-and-Run)\n";
        std::cout << "  7. Barrier Methods (Dikin, John, Vaidya, CRHMC)\n";
        std::cout << "  8. Boundary & Other Methods (Shake & Bake)\n";
        std::cout << "  0. Return to Main Menu\n";

        std::cout << "\n" << std::string(60, '=') << "\n";

        // Display the selected information block
        if (choice == -1) {
            std::cout << "Please select an option from the menu above.\n";
        } else {
            switch (choice) {
                case 1:
                    std::cout << R"(--- Markov Chain Monte Carlo (MCMC) ---
MCMC is a class of algorithms used to sample from a probability 
distribution. Because high-dimensional polytopes are difficult 
to sample from directly, MCMC constructs a random walk (a 
Markov chain) where each step depends only on the previous one. 
Over time, the steps of this walk distribute themselves 
according to the target distribution (usually uniform for 
polytopes).)" << '\n';
                    break;
                case 2:
                    std::cout << R"(--- Effective Sample Size (ESS) ---
Since MCMC samples are generated via a random walk, consecutive 
samples are correlated (they are near each other in space). 
ESS calculates how many statistically independent samples your 
autocorrelated chain is equivalent to. 

For example, 10,000 raw samples might only yield an ESS of 500, 
meaning the walk is moving slowly through the space. A higher 
ESS indicates better mixing and a more efficient algorithm.)" << '\n';
                    break;
                case 3:
                    std::cout << R"(--- Potential Scale Reduction Factor (PSRF) ---
Also known as the Gelman-Rubin diagnostic. It is used to evaluate 
if your MCMC walk has converged to the target distribution. 

It compares the variance between multiple independent random 
walks to the variance within those individual walks. A PSRF 
value close to 1.0 (typically < 1.1) suggests that the chains 
have converged and are exploring the same space properly.)" << '\n';
                    break;
                case 4:
                    std::cout << R"(--- Kolmogorov-Smirnov (KS) Test ---
The KS test is a statistical test used to compare a sample 
distribution with a reference probability distribution, or to 
compare two sample distributions. 

In this benchmark, it is often used to check the 1-dimensional 
marginals of the sampled points to see if they match the 
theoretically expected uniform distribution of the polytope, 
verifying the correctness of the walk.)" << '\n';
                    break;
                case 5:
                    std::cout << R"(--- Polytope Rounding ---
Polytopes can be highly skewed, elongated, or "skinny". Random 
walks struggle in skinny spaces because they easily get stuck 
bouncing off the walls, leading to high autocorrelation.

Rounding computes an ellipsoid (e.g., Maximum Volume Inscribed 
Ellipsoid) that fits inside the polytope, and uses it to apply 
a linear transformation. This transformation makes the polytope 
look more like a sphere (isotropic). The random walk runs much 
faster in this rounded space, and the generated samples are then 
transformed back to the original space.)" << '\n';
                    break;
                case 6:
                    std::cout << R"(--- Geometric Methods ---
These methods rely purely on the geometric properties of the space:
- Hit-and-Run: Picks a random direction, computes the line segment 
  inside the polytope along that direction, and jumps to a uniform 
  random point on that segment.
- Ball Walk: Proposes a point uniformly within a small ball around 
  the current point. Rejects the step if it falls outside the polytope.
- Billiards: Simulates particle trajectories bouncing off the 
  polytope's inner boundaries, like a billiard ball.)" << '\n';
                    break;
                case 7:
                    std::cout << R"(--- Barrier Methods ---
These methods utilize barrier functions (often borrowed from 
Interior Point methods in optimization) to keep the walk away from 
the boundaries:
- Dikin Walk: Uses the log-barrier Hessian to define an ellipsoid 
  around the current point that is guaranteed to be inside the 
  polytope, proposing the next step from within this ellipsoid.
- Vaidya / John Walks: Use more complex barrier functions 
  (Vaidya's barrier or John's ellipsoid) to take larger steps 
  without hitting the boundary.
- CRHMC: Simulates a particle moving inside the polytope under 
  the Hamiltonian laws of motion, where the boundaries act as 
  immovable walls.)" << '\n';
                    break;
                case 8:
                    std::cout << R"(--- Boundary & Other Methods ---
While standard methods sample the interior volume of the polytope, 
some tasks require sampling exactly from the surface (the facets).

- Shake & Bake: An algorithm specifically designed to generate 
  points uniformly distributed on the boundary of a polytope. It 
  "shakes" to find a direction and "bakes" by moving along the 
  boundary facets.)" << '\n';
                    break;
            }
        }
        std::cout << std::string(60, '=') << "\n";

        // Prompt the user for the next action
        choice = get_valid_int("Enter your choice (0-8): ", 0, 8);

        if (choice == 0) {
            // Clear the screen one last time before returning to the main menu
#if defined(_WIN32)
            int final_ret = std::system("cls");
#else
            int final_ret = std::system("clear");
#endif
            (void)final_ret; // Silence the warning
            return;
        }
    }
}

// Handles the setup flow when the user chooses Option 1
inline void setup_benchmark_options(BenchmarkConfig& config, std::string& walk_choice) {
    std::cout << "\n--- Benchmark Setup (Type 'q' to quit at any prompt) ---\n";

    // 1. Polytope Choice
    std::string poly_prompt = "Select a Polytope:\n"
                              "1. Cube\n"
                              "2. Simplex\n"
                              "3. Birkhoff\n"
                              "4. Cross\n"
                              "5. OrderPolytope\n"
                              "6. Custom\n"
                              "Choice (1-6): ";
    
    int p_choice = get_valid_int(poly_prompt, 1, 6);
    switch (p_choice) {
        case 1: config.polytope_choice = "Cube"; break;
        case 2: config.polytope_choice = "Simplex"; break;
        case 3: config.polytope_choice = "Birkhoff"; break;
        case 4: config.polytope_choice = "Cross"; break;
        case 5: config.polytope_choice = "OrderPolytope"; break;
        case 6: 
            config.polytope_choice = "Custom"; 
            std::cout << "\n[Note] You selected Custom. The tool will use the 'custom_A_file'\n"
                      << "and 'custom_b_file' exactly as defined in your JSON config.\n";
            break;
    }

    // 2. Dimension (Skip if Custom Polytope)
    if (config.polytope_choice != "Custom") {
        std::cout << "\n";
        // Directly overwrite the JSON dimensions array with the user's list
        config.dimensions = get_valid_dimension_list("Enter dimensions separated by commas (e.g., 10,20,50): ");
    } else {
        config.dimensions.clear();
        config.dimensions.push_back(0); 
    }

    // 3. Stopping Criterion
    std::cout << "\n";
    std::string stop_prompt = "Select stopping criterion:\n"
                              "1. Target ESS\n"
                              "2. Fixed number of samples\n"
                              "Choice (1-2): ";
    int stop_choice = get_valid_int(stop_prompt, 1, 2);

    if (stop_choice == 1) {
        config.use_dynamic_batch = true;
        config.target_ESS = get_valid_int("\nEnter Target ESS (10 to 100000): ", 10, 100000);
    } else {
        config.use_dynamic_batch = false;
        int fixed_samples = get_valid_int("\nEnter number of samples to generate (1 to 100000): ", 1, 100000);
        
        // Override the 'samples' value for all walks in the configuration.
        // This ensures your get_initial_batch_size() function reads the correct user input.
        for (auto& pair : config.walk_settings) {
            pair.second.samples = fixed_samples;
        }
    }

    // 4. Time Limit
    std::cout << "\n";
    int time_lim = get_valid_int("Enter time limit in seconds (1 to 86400): ", 1, 86400);
    config.time_limit_sec = static_cast<double>(time_lim);

    // 5. Base Seed
    std::cout << "\n";
    config.base_seed = get_valid_int("Enter base seed (e.g., 42): ", 0, 2147483647);

    // 6. Rotation Angle
    std::cout << "\n";
    config.angle = get_valid_int("Enter rotation angle in degrees (0 to 360): ", 0, 360);

    // 7. Write to File
    std::cout << "\n";
    int write_choice = get_valid_int("Write samples to file?\n1. Yes\n2. No\nChoice (1-2): ", 1, 2);
    config.write_to_file = (write_choice == 1);

    // 8. Rounding
    std::cout << "\n";
    std::cout << "Select a rounding method:\n";
    std::cout << "1. None (Disabled)\n";
    std::cout << "2. Max Ellipsoid\n";
    std::cout << "3. Log Barrier\n";
    std::cout << "4. Vaidya Barrier\n";
    std::cout << "5. Volumetric Barrier\n";
    int r_choice = get_valid_int("Choice (1-5): ", 1, 5);

    switch (r_choice) {
        case 1: 
            config.rounding = false; 
            config.rounding_method = "none";
            break;
        case 2: 
            config.rounding = true; 
            config.rounding_method = "max_ellipsoid"; 
            break;
        case 3: 
            config.rounding = true; 
            config.rounding_method = "log_barrier"; 
            break;
        case 4: 
            config.rounding = true; 
            config.rounding_method = "vaidya_barrier"; 
            break;
        case 5: 
            config.rounding = true; 
            config.rounding_method = "volumetric_barrier"; 
            break;
    }
    // 9. Method Choice
    std::cout << "\n";
    std::string strategy_prompt = "Select a Walk Method Strategy:\n"
                                  "1. All (Run all enabled in JSON)\n"
                                  "2. Auto (Select based on dimension)\n"
                                  "3. Manual Selection\n"
                                  "Choice (1-3): ";
    
    int strat_choice = get_valid_int(strategy_prompt, 1, 3);

    if (strat_choice == 1) {
        walk_choice = "All"; 
        config.auto_walk = false;
    } 
    else if (strat_choice == 2) {
        walk_choice = "All"; 
        config.auto_walk = true; 
    } 
    else {
        config.auto_walk = false;
        std::string dist_prompt = "\nSelect Target Distribution:\n"
                                  "1. Uniform\n"
                                  "2. Exponential (Gaussian)\n"
                                  "Choice (1-2): ";
        int dist_choice = get_valid_int(dist_prompt, 1, 2);

        if (dist_choice == 2) {
            std::string exp_prompt = "\nSelect Exponential (Gaussian) Method:\n"
                                     "1. Gaussian Ball\n"
                                     "2. Gaussian Coordinate Direction Hit and Run\n"
                                     "Choice (1-2): ";
            int exp_choice = get_valid_int(exp_prompt, 1, 2);
            switch(exp_choice) {
                case 1: walk_choice = "GaussianBallWalk"; break;
                case 2: walk_choice = "GaussianCDHRWalk"; break;
            }
        } 
        else {
            std::string cat_prompt = "\nSelect Uniform Method Category:\n"
                                     "1. Geometric (Ball, Billiards, Hit-and-Run)\n"
                                     "2. Barrier (Dikin, John, Vaidya, CRHMC)\n"
                                     "3. Boundary & Other (Shake & Bake)\n"
                                     "Choice (1-3): ";
            int cat_choice = get_valid_int(cat_prompt, 1, 3);

            if (cat_choice == 1) {
                std::string geom_prompt = "\nSelect Geometric Method:\n"
                                          "1. Ball\n"
                                          "2. Billiard\n"
                                          "3. Accelerated Billiard\n"
                                          "4. Sparse Billiard\n"
                                          "5. Coordinate Direction Hit and Run (CDHR)\n"
                                          "6. Random Direction Hit and Run (RDHR)\n"
                                          "Choice (1-6): ";
                int geom_choice = get_valid_int(geom_prompt, 1, 6);
                switch(geom_choice) {
                    case 1: walk_choice = "BallWalk"; break;
                    case 2: walk_choice = "BilliardWalk"; break;
                    case 3: walk_choice = "AcceleratedBilliardWalk"; break;
                    case 4: walk_choice = "SparseBilliardWalk"; break;
                    case 5: walk_choice = "CDHRWalk"; break;
                    case 6: walk_choice = "RDHRWalk"; break;
                }
            } 
            else if (cat_choice == 2) {
                std::string bar_prompt = "\nSelect Barrier Method:\n"
                                         "1. Dikin\n"
                                         "2. John\n"
                                         "3. Vaidya\n"
                                         "4. Constrained Riemannian Hamiltonian Monte Carlo (CRHMC)\n"
                                         "Choice (1-4): ";
                int bar_choice = get_valid_int(bar_prompt, 1, 4);
                switch(bar_choice) {
                    case 1: walk_choice = "DikinWalk"; break;
                    case 2: walk_choice = "JohnWalk"; break;
                    case 3: walk_choice = "VaidyaWalk"; break;
                    case 4: walk_choice = "CRHMCWalk"; break;
                }
            } 
            else {
                std::string other_prompt = "\nSelect Boundary & Other Method:\n"
                                           "1. Billiard Shake and Bake\n"
                                           "2. Shake and Bake\n"
                                           "3. Boundary Coordinate Direction Hit and Run\n"
                                           "4. Boundary Random Direction Hit and Run\n"
                                           "Choice (1-4): ";
                int other_choice = get_valid_int(other_prompt, 1, 4);
                switch(other_choice) {
                    case 1: walk_choice = "BilliardShakeAndBakeWalk"; break;
                    case 2: walk_choice = "ShakeAndBakeWalk"; break;
                    case 3: walk_choice = "BCDHRWalk"; break;
                    case 4: walk_choice = "BRDHRWalk"; break;
                }
            }
        }
    }

    std::cout << "\nSetup complete! Proceeding to benchmark...\n";
}

// Returns true if the benchmark should proceed, false if the user chose to quit.
inline bool run_interactive_menu(BenchmarkConfig& config, std::string& walk_choice) {
    std::cout << "\n==========================================\n";
    std::cout << "   Welcome to the Volesti Benchmark Tool  \n";
    std::cout << "==========================================\n\n";

    while (true) {
        std::cout << "Main Menu:\n";
        std::cout << "1. Start benchmark\n";
        std::cout << "2. Show help\n";
        std::cout << "3. Quit\n";
        
        int choice = get_valid_int("Enter your choice (1-3): ", 1, 3);

        switch (choice) {
            case 1:
                setup_benchmark_options(config, walk_choice);
                return true; 
            case 2:
                show_help();
                break; 
            case 3:
                std::cout << "\nExiting program. Goodbye!\n";
                return false; 
        }
    }
}