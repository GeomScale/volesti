// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2026 Vissarion Fisikopoulos
// Copyright (c) 2018-2026 Apostolos Chalkis
// Copyright (c) 2026      Dimitrios Pavlou

// Contributed and/or modified by Dimitrios Pavlou, as part of Google Summer of Code 2026 program

// Licensed under GNU LGPL.3, see LICENCE file

#include <vector>
#include "doctest.h"
#include "Highs.h"

void call_test_small_lp()
{
    // Small LP: max x+y s.t. x+y <= 2, x,y>= 0
    // Optimal: x=1, y=1, obj=2
    Highs highs;
    highs.setOptionValue("output_flag", false);

    highs.changeObjectiveSense(ObjSense::kMaximize);

    highs.addVar(0.0, kHighsInf);
    highs.addVar(0.0, kHighsInf);
    highs.changeColCost(0, 1);
    highs.changeColCost(1, 1);

    std::vector<HighsInt>indices = {0, 1};
    std::vector<double> values = {1, 1};
    highs.addRow(-kHighsInf, 2.0, 2, indices.data(), values.data());
    
    highs.run();

    CHECK(highs.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(std::abs(highs.getObjectiveValue()-2.0) < 1e-10);
}

void call_test_infeasible_lp() 
{
    // LP: min 0 s.t. x >= 1 x <= -1
    // Infeasible LP
    Highs highs;
    highs.setOptionValue("output_flag", false);

    highs.addVar(0.0, kHighsInf);
    std::vector<HighsInt> indices = {0};
    std::vector<double> values = {1.0};

    highs.addRow(1.0, kHighsInf, 1, indices.data(), values.data());
    highs.addRow(-kHighsInf, -1.0, 1, indices.data(), values.data());

    highs.run();

    CHECK(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}

void call_test_unbounded_lp() 
{
    // LP: min -x s.t. x >= 0
    // Unbounded LP
    Highs highs;
    highs.setOptionValue("output_flag", false);

    highs.addVar(0.0, kHighsInf);
    highs.changeColCost(0, -1.0);

    highs.run();

    CHECK(highs.getModelStatus() == HighsModelStatus::kUnbounded);
}

void call_test_max_inner_ball_lp(unsigned d) 
{
    // LP: max r s.t. a^Tx+r||ai|| <= bi
    // For the hypercube the constraints are (+/-)x+r<=1
    // Optimal: r=1 x=(0,0,...,0)
    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.changeObjectiveSense(ObjSense::kMaximize);

    for (unsigned i = 0; i < d; ++i)
        highs.addVar(-kHighsInf, kHighsInf);

    highs.addVar(0.0, kHighsInf);
    highs.changeColCost(d, 1.0);

    for (unsigned i = 0; i < d; ++i) {
        std::vector<HighsInt> indices = {HighsInt(i), HighsInt(d)};
        std::vector<double> values = {1.0, 1.0};
        highs.addRow(-kHighsInf, 1.0, 2, indices.data(), values.data());

        values[0] *= -1;
        highs.addRow(-kHighsInf, 1.0, 2, indices.data(), values.data());
    }

    highs.run();

    CHECK(highs.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(std::abs(highs.getObjectiveValue() -1.0) < 1e-10);

    const HighsSolution& sol = highs.getSolution();
    for (unsigned i = 0; i < d; ++i)
        CHECK(std::abs(sol.col_value[i]) < 1e-10);

    CHECK(std::abs(sol.col_value[d]-1.0) < 1e-10);
}

TEST_CASE("test_small_lp") {
    call_test_small_lp();
}

TEST_CASE("test_infeasible_lp") {
    call_test_infeasible_lp();
}

TEST_CASE("test_unbounded_lp") {
    call_test_unbounded_lp();
}

TEST_CASE("test_max_inner_ball_lp") {
    call_test_max_inner_ball_lp(3);
    call_test_max_inner_ball_lp(5);
    call_test_max_inner_ball_lp(10);
    call_test_max_inner_ball_lp(20);
}