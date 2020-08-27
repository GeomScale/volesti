// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "NonSymmetricIPM.h"
#include "MonomialsClass.h"
#include "EnvelopeProblemSDP.h"
#include "EnvelopeProblemSOS.h"
#include "spdlog/spdlog.h"
#include "spdlog/cfg/env.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <fstream>
#include "misc/tests.cpp"

int main(int const argc, char **argv) {

    srand(time(nullptr));

    auto console = spdlog::stdout_color_mt("console");
    console->info("Logger level is {}", console->level());
//    Eigen::setNbThreads(1);
    console->info("Num threads: {}", Eigen::nbThreads());

    std::ifstream instance_file;
    std::ifstream config_file;

    std::string instance_file_str;
    if (argc < 2) {
        console->info("No data file provided. The default file will be used instead.");
        instance_file_str = "../config/instance.json";
        instance_file.open(instance_file_str);
        if (not instance_file.is_open()) {
            instance_file_str = "config/instance.json";
            instance_file.open(instance_file_str);
        }
        if (not instance_file.is_open()) {
            console->error("Could not locate file.");
            return 1;
        }
    } else {
        instance_file.open(argv[1]);
        if (not instance_file.is_open()) {
            console->error("Could not locate file {}", argv[1]);
            return 1;
        }
    }

    std::string config_file_str;
    if (argc < 3) {
        console->info("No configuration file provided. The default file will be used instead.");
        config_file_str = "../config/config.json";
        config_file.open(config_file_str);
        if (not config_file.is_open()) {
            config_file_str = "config/config.json";
            config_file.open(config_file_str);
        }
        if (not config_file.is_open()) {
            console->error("Could not locate file.");
            return 1;
        }
    } else {
        config_file_str = argv[2];
        config_file.open(config_file_str);
        if (not config_file.is_open()) {
            console->error("Could not locate file {}", argv[2]);
            return 1;
        }
    }

    bool run_tests = false;
    for (unsigned i = 2; i < (unsigned) argc; i++) {
        std::string arg_str(argv[i]);
        if (arg_str == "run_tests") {
            run_tests = true;
        }
    }

    if (run_tests) {
        console->info("Test LP and SDP solver");
        assert(test_lp_solver_random(2, 5));
        assert(test_sdp_solver_random_lp_formulation(2, 5));
        console->info("Tests completed.");
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        test_sdp_solver(config_file_str);
        //Reset file to be read again by proper SOS solver.
        config_file.clear();
        config_file.seekg(0);
    }

    EnvelopeProblemSOS envelopeProblemSos(instance_file_str, config_file_str);
    Instance instance_interp = envelopeProblemSos.construct_SOS_instance();
    NonSymmetricIPM sos_solver_interp(instance_interp, config_file_str);
    sos_solver_interp.run_solver();
    envelopeProblemSos.print_solution(sos_solver_interp.get_solution());
    envelopeProblemSos.plot_polynomials_and_solution(sos_solver_interp.get_solution());
    return 0;
}

