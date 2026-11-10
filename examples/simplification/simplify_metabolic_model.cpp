// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <filesystem>
#include <string>
#include <ctime>
#include <iomanip>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "preprocess/constrained_polytope/bigg_parser.hpp"
#include "preprocess/constrained_polytope/simplifier.hpp"
#include "preprocess/constrained_polytope/transformation.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef ConstrainedPolytope<Point> Polytope;

// Prints a label and its value (a row of the output).
// @tparam T the value type
// @param label the label
// @param value the value
template <typename T>
void row(std::string const& label, T const& value) {
    std::cout << "  " << std::left << std::setw(20) << label << value << "\n"; 
}

// Prints the statistics of a simplification run.
// @param name the name of the model
// @param P the original polytope
// @param res the simplification result
// @param trans_t the cpu time of the transformation in seconds
void print_statistics(std::string const& name,
                      Polytope const& P,
                      SimplifierResult<Point> const& res,
                      double trans_t
                     )
{
    std::string const line(48, '=');
    Polytope const& Ps = res.polytope;

    std::cout << line << "\n" << name << "\n" << line << "\n";

    std::cout << "\nPolytope\n";
    row("variables", P.getDimension());
    row("equalities", std::to_string(P.getNumEqualities())+" -> "+std::to_string(Ps.getNumEqualities()));
    row("finite bounds", std::to_string(P.getNumFiniteBounds())+" -> "+std::to_string(Ps.getNumFiniteBounds()));
    unsigned const before = P.getNumFiniteBounds();
    unsigned const after = Ps.getNumFiniteBounds();
    row("removed", std::to_string(before ? 100.0*(before-after)/before : 0.0)+"%");

    std::cout << "\nStages\n";
    row("pinned variables", res.presolve.pinned_vars);
    row("fixed variables", res.reduce.fixed_vars);
    row("reduce passes", res.reduce.passes);
    row("removed rows", res.dependent_rows.removed_rows);
    row("interior slack", res.interior_point.slack);
    row("total LPs", res.solved_lps);

    std::cout << "\nTime (cpu)\n";
    row("simplification", std::to_string(res.times.total())+" s");
    row("transformation", std::to_string(trans_t)+" s");

    std::cout << "\n" << line << "\n status: " << status_to_string(res.status) << "\n" << line << "\n";
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " <model.json>" << std::endl;
        return 1;
    }

    std::filesystem::path const input = argv[1];
    if (!std::filesystem::exists(input)) {
        std::cerr << "input not found: " << input << std::endl;
        return 1;
    }

    std::string const name = input.stem().string();

    Polytope P = parse_from_json<Point>(input.string());

    SimplifierConfig config;
    auto res = simplify(P, config);

    if (res.status != SimplifierStatus::OK) {
        print_statistics(name, P, res, 0.0);
        return 2;
    }

    // Grabs the full dimensional polytope
    std::clock_t c0 = std::clock();
    auto [HP, shift, N] = transform(res.polytope);
    double const trans_t = double(std::clock()-c0)/CLOCKS_PER_SEC;

    print_statistics(name, P, res, trans_t);

    // The full dimensional polytope is {x : A x <= b}, mapped back by y = N x + shift
    std::cout << "full dimensional polytope: dimension " << HP.dimension()
              << ", " << HP.num_of_hyperplanes() << " facets"
              << " (N is " << N.rows() << "x" << N.cols()
              << ", shift has " << shift.size() << " entries)" << std::endl;
    return 0;
}