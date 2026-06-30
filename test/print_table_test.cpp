// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis
// Copyright (c) 2021 Maios Papachristou

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "misc/print_table.hpp"
#include <sstream>
#include <string>

TEST_CASE("VariadicTable::empty_table")
{
    VariadicTable<std::string, int> vt({"Name", "Age"});
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    CHECK(output.size() > 0);
    // Should contain the headers
    CHECK(output.find("Name") != std::string::npos);
    CHECK(output.find("Age") != std::string::npos);
}

TEST_CASE("VariadicTable::single_row")
{
    VariadicTable<std::string, int> vt({"Name", "Age"});
    vt.addRow("Alice", 30);
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    CHECK(output.size() > 0);
    CHECK(output.find("Alice") != std::string::npos);
    CHECK(output.find("30") != std::string::npos);
}

TEST_CASE("VariadicTable::multiple_rows_and_types")
{
    VariadicTable<std::string, double, int> vt({"Item", "Price", "Qty"});
    vt.addRow("Apple", 1.25, 10);
    vt.addRow("Banana", 0.75, 20);
    vt.addRow("Cherry", 2.50, 15);
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    CHECK(output.find("Apple") != std::string::npos);
    CHECK(output.find("Banana") != std::string::npos);
    CHECK(output.find("Cherry") != std::string::npos);
    CHECK(output.find("1.25") != std::string::npos);
    CHECK(output.find("0.75") != std::string::npos);
    CHECK(output.find("2.5") != std::string::npos);
    CHECK(output.find("10") != std::string::npos);
    CHECK(output.find("20") != std::string::npos);
    CHECK(output.find("15") != std::string::npos);
}

TEST_CASE("VariadicTable::column_format_scientific")
{
    VariadicTable<std::string, double> vt({"Label", "Value"});
    vt.addRow("Pi", 3.14159);
    vt.addRow("E", 2.71828);
    vt.setColumnFormat({VariadicTableColumnFormat::AUTO,
                        VariadicTableColumnFormat::SCIENTIFIC});
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    // Scientific format should contain 'e'
    bool has_exp = (output.find("e+") != std::string::npos) || (output.find("e-") != std::string::npos);
    CHECK(has_exp);
}

TEST_CASE("VariadicTable::column_format_fixed")
{
    VariadicTable<std::string, double> vt({"Label", "Value"});
    vt.addRow("Half", 0.5);
    vt.addRow("Third", 1.0 / 3.0);
    vt.setColumnFormat({VariadicTableColumnFormat::AUTO,
                        VariadicTableColumnFormat::FIXED});
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    // FIXED keeps trailing zeros via precision, so we should find the decimal point
    CHECK(output.find(".") != std::string::npos);
    // Both values should appear
    CHECK(output.find("0.5") != std::string::npos);
}

TEST_CASE("VariadicTable::column_format_percent")
{
    VariadicTable<std::string, double> vt({"Label", "Pct"});
    vt.addRow("A", 0.1234);
    vt.addRow("B", 1.0);
    vt.setColumnFormat({VariadicTableColumnFormat::AUTO,
                        VariadicTableColumnFormat::PERCENT});
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    // PERCENT forces setprecision(2) and width 6, so 0.1234 becomes "0.12" and 1.0 becomes "1.00"
    CHECK(output.find("1.00") != std::string::npos);
    CHECK(output.find("0.12") != std::string::npos);
}

TEST_CASE("VariadicTable::column_precision")
{
    VariadicTable<std::string, double> vt({"Label", "Value"});
    vt.addRow("Pi", 3.1415926535);
    vt.setColumnPrecision({0, 3});
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    // With defaultfloat, setprecision(3) means 3 significant digits, so 3.14159 -> "3.14"
    CHECK(output.find("3.14") != std::string::npos);
}

TEST_CASE("VariadicTable::static_column_size")
{
    // Without static_column_size a plain type like double that is not
    // integral and has no .size() gets width 0, so set it explicitly.
    VariadicTable<double, int> vt({"X", "Y"}, 10);
    vt.addRow(1.234, 42);
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    CHECK(output.find("1.234") != std::string::npos);
    CHECK(output.find("42") != std::string::npos);
}

TEST_CASE("VariadicTable::custom_cell_padding")
{
    VariadicTable<std::string, int> vt({"Name", "Age"}, 0, 3);
    vt.addRow("Bob", 25);
    std::ostringstream oss;
    vt.print(oss);
    std::string output = oss.str();
    // With padding 3 there should be at least 3 spaces between '|' and text
    // (the value is right-justified for arithmetic types, left-justified for strings)
    bool bob_found = output.find("Bob") != std::string::npos;
    bool twentyfive_found = output.find("25") != std::string::npos;
    CHECK(bob_found);
    CHECK(twentyfive_found);
}
