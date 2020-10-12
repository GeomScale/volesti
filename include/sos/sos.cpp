// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "NonSymmetricIPM.h"
#include "../../examples/EnvelopeProblemSOS/EnvelopeProblemSOS.h"
#include "spdlog/spdlog.h"
#include "spdlog/cfg/env.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <fstream>

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
        instance_file_str = argv[1];
        instance_file.open(std::string(instance_file_str));
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

    EnvelopeProblemSOS<double> envelopeProblemSos(instance_file_str, config_file_str);
    Instance<double> instance_interp = envelopeProblemSos.construct_SOS_instance();
    NonSymmetricIPM<double> sos_solver_interp(instance_interp, config_file_str);
    int outcome = sos_solver_interp.run_solver();

    if(outcome == NonSymmetricIPM<double>::Termination::SUCCESS){
        envelopeProblemSos.print_solution(sos_solver_interp.get_solution());
        envelopeProblemSos.plot_polynomials_and_solution(sos_solver_interp.get_solution());
    } else if (sos_solver_interp._type_cast_if_unsuccessful) {
        console->error("Switch to more precise double type");
        NonSymmetricIPM<long double> *long_double_solver = sos_solver_interp.cast_with_product_barrier<long double>();
        long_double_solver->run_solver();
        Solution<long double> long_double_sol = long_double_solver->get_solution();
        envelopeProblemSos.print_solution(long_double_sol.cast<double>());
        envelopeProblemSos.plot_polynomials_and_solution(long_double_sol.cast<double>());
    }

    return 0;
}

