#include "../include/walk_parameters.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>

BenchmarkConfig load_benchmark_config(const std::string& filepath) {
    BenchmarkConfig config;
    boost::property_tree::ptree pt;

    try {
        // Parse the JSON file
        boost::property_tree::read_json(filepath, pt);

        // Read Global Settings
        config.target_ESS      = pt.get<unsigned int>("global_settings.target_ESS", 500);
        config.time_limit_sec  = pt.get<double>("global_settings.time_limit_sec", 1200.0);
        config.base_seed       = pt.get<int>("global_settings.base_seed", 42);

        //dimensions array 
        auto dims_node = pt.get_child_optional("global_settings.dimensions");
        if (dims_node) {
            for (const auto& item : *dims_node) {
                config.dimensions.push_back(item.second.get_value<unsigned int>());
            }
        } else {
            // fallback for older configs
            config.dimensions.push_back(pt.get<unsigned int>("global_settings.dimensions", 100));
        }

        config.polytope_choice = pt.get<std::string>("global_settings.polytope_choice", "Cube");
        config.custom_A_file   = pt.get<std::string>("global_settings.custom_A_file", "");
        config.custom_b_file   = pt.get<std::string>("global_settings.custom_b_file", "");
        config.angle           = pt.get<unsigned int>("global_settings.rotation_angle", 0);
        config.use_dynamic_batch = pt.get<bool>("global_settings.dynamic_batch_size", true);
        config.write_to_file     = pt.get<bool>("global_settings.write_to_file", false);
        config.rounding          = pt.get<bool>("global_settings.rounding", false);
        config.rounding_method   = pt.get<std::string>("global_settings.rounding_method", "max_ellipsoid");
        config.auto_walk         = pt.get<bool>("global_settings.auto_walk", false);
        config.show_console_logs = pt.get<bool>("global_settings.show_console_logs", true);
        config.show_menu         = pt.get<bool>("global_settings.show_menu", false);

        // Iterate through the walks JSON object
        for (const auto& walk_node : pt.get_child("walks")) {
            std::string walk_name = walk_node.first;
            WalkSettings settings;

            settings.enabled             = walk_node.second.get<bool>("enabled", true);
            settings.samples             = walk_node.second.get<unsigned int>("samples", 1000);
            settings.walk_len_multiplier = walk_node.second.get<unsigned int>("walk_len_multiplier", 0);
            settings.walk_len_base       = walk_node.second.get<unsigned int>("walk_len_base", 1);

            //For Gaussian walks
            settings.a_i_param           = walk_node.second.get<double>("a_i_param", 1.0);

            config.walk_settings[walk_name] = settings;
        }
    } catch (const boost::property_tree::ptree_error& e) {
        std::cerr << "Error reading config file: " << e.what() << "\n";
        std::cerr << "Falling back to hardcoded defaults.\n";
    }

    return config;
}

unsigned int get_initial_batch_size(const std::string& walk_name, const BenchmarkConfig& config) {
    auto it = config.walk_settings.find(walk_name);
    if (it != config.walk_settings.end()) {
        return it->second.samples + (2 * config.dimension);
    }
    return 1000; 
}

unsigned int compute_dynamic_walk_len(const std::string& walk_name, unsigned int dim, const BenchmarkConfig& config) {
    auto it = config.walk_settings.find(walk_name);
    if (it != config.walk_settings.end()) {
        return (dim * it->second.walk_len_multiplier) + it->second.walk_len_base;
    }
    return 1; 
}