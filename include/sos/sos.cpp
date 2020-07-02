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
#include "tests.cpp"

int main(int const argc, char **argv) {

    srand(time(nullptr));

    auto console = spdlog::stdout_color_mt("console");
    console->info("Logger level is {}", console->level());

    std::ifstream file;
    if (argc < 2) {
        console->info("No data file provided. The default file will be used instead.");
        file.open("../config/default.txt");
        if (not file.is_open()) {
            file.open("config/default.txt");
        }
        if (not file.is_open()) {
            console->error("Could not locate file.");
            return 1;
        }
    } else {
        file.open(argv[1]);
        if(not file.is_open()){
            console->error("Could not local file {}", argv[1]);
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
    }

    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    int max_degree;
    iss >> max_degree;

    HyperRectangle hyperRectangle;
    //Note: Keep interval bounds for now.
    IPMDouble const interval_lower_bound = -1.;
    IPMDouble const interval_upper_bound = 1.;
    hyperRectangle.push_back(std::pair<IPMDouble, IPMDouble>(interval_lower_bound,
                                                             interval_upper_bound));
    EnvelopeProblemSOS envelopeProblemSos(1, max_degree, hyperRectangle);

    while (std::getline(file, line)) {
        std::istringstream poly_stream(line);
        InterpolantVector sos_poly = envelopeProblemSos.generate_zero_polynomial();
        IPMDouble val;
        unsigned idx = 0;
        while (poly_stream >> val) {
            if (idx >= sos_poly.size()) {
                return 1;
            }
            sos_poly[idx++] = val;
        }
        envelopeProblemSos.add_polynomial(sos_poly);
        console->info("Polynomial added.");
    }

    Instance instance_interp = envelopeProblemSos.construct_SOS_instance();

    NonSymmetricIPM sos_solver_interp(instance_interp);

    sos_solver_interp.run_solver();

    envelopeProblemSos.print_solution(sos_solver_interp.get_solution());
    envelopeProblemSos.plot_polynomials_and_solution(sos_solver_interp.get_solution());
    return 0;
}

