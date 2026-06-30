// VolEsti (volume computation and sampling library)
// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include "misc/poset.h"
#include "misc/misc.h"
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

// =============================================================================
// Poset tests
// =============================================================================

TEST_CASE("Poset default constructor") {
    Poset p;
    CHECK(p.num_elem() == 0);
    CHECK(p.num_relations() == 0);
}

TEST_CASE("Poset verify valid") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 1}, {1, 2}
    };
    auto result = Poset::verify(rel, 3);
    CHECK(result.size() == 2);
    CHECK(result[0].first == 0);
    CHECK(result[0].second == 1);
}

TEST_CASE("Poset verify empty relations") {
    std::vector<std::pair<unsigned int, unsigned int>> rel;
    auto result = Poset::verify(rel, 5);
    CHECK(result.empty());
}

TEST_CASE("Poset verify single element") {
    std::vector<std::pair<unsigned int, unsigned int>> rel;
    auto result = Poset::verify(rel, 1);
    CHECK(result.empty());
}

TEST_CASE("Poset verify out of range element") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 5}  // 5 >= n (3)
    };
    CHECK_THROWS(Poset::verify(rel, 3));
}

TEST_CASE("Poset verify out of range negative-like") {
    // Test with a large value beyond n
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {2, 3}  // 3 is out of range when n=3
    };
    CHECK_THROWS(Poset::verify(rel, 3));
}

TEST_CASE("Poset verify cycle") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 1}, {1, 2}, {2, 0}  // cycle 0->1->2->0
    };
    CHECK_THROWS(Poset::verify(rel, 3));
}

TEST_CASE("Poset verify self-loop") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 0}  // self-loop creates a cycle
    };
    CHECK_THROWS(Poset::verify(rel, 3));
}

TEST_CASE("Poset valid constructor") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 1}, {0, 2}, {1, 3}, {2, 3}
    };
    Poset p(4, rel);
    CHECK(p.num_elem() == 4);
    CHECK(p.num_relations() == 4);
}

TEST_CASE("Poset invalid constructor throws") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 1}, {1, 0}  // cycle
    };
    CHECK_THROWS(Poset(2, rel));
}

TEST_CASE("Poset get_relation") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}, {2, 3}};
    Poset p(4, rel);
    auto r0 = p.get_relation(0);
    CHECK(r0.first == 0);
    CHECK(r0.second == 1);
    auto r1 = p.get_relation(1);
    CHECK(r1.first == 2);
    CHECK(r1.second == 3);
}

TEST_CASE("Poset num_elem and num_relations") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}};
    Poset p(3, rel);
    CHECK(p.num_elem() == 3);
    CHECK(p.num_relations() == 1);
}

TEST_CASE("Poset is_in with double") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}, {1, 2}};
    Poset p(3, rel);
    // x0 <= x1 <= x2: e.g., [0.1, 0.2, 0.3] satisfies
    std::vector<double> pt1 = {0.1, 0.2, 0.3};
    CHECK(p.is_in(pt1));
    // x0 > x1 violates: [0.5, 0.2, 0.3]
    std::vector<double> pt2 = {0.5, 0.2, 0.3};
    CHECK(!p.is_in(pt2));
}

TEST_CASE("Poset is_in with tolerance") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}};
    Poset p(2, rel);
    // Nearly equal with tolerance
    std::vector<double> pt = {1.0001, 1.0};
    CHECK(!p.is_in(pt));          // strict
    CHECK(p.is_in(pt, 0.001));    // with tolerance
}

TEST_CASE("Poset is_in empty poset always true") {
    std::vector<std::pair<unsigned int, unsigned int>> empty_rel;
    Poset p(5, empty_rel);
    std::vector<double> pt = {100.0, -200.0, 0.0, 1.5, -3.2};
    CHECK(p.is_in(pt));
}

TEST_CASE("Poset topologically_sorted_list") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {
        {0, 1}, {0, 2}, {1, 3}, {2, 3}
    };
    Poset p(4, rel);
    auto sorted = p.topologically_sorted_list();
    CHECK(sorted.size() == 4);
    // 0 must come before 1 and 2
    // 1 and 2 must come before 3
    auto idx0 = std::find(sorted.begin(), sorted.end(), 0u) - sorted.begin();
    auto idx1 = std::find(sorted.begin(), sorted.end(), 1u) - sorted.begin();
    auto idx2 = std::find(sorted.begin(), sorted.end(), 2u) - sorted.begin();
    auto idx3 = std::find(sorted.begin(), sorted.end(), 3u) - sorted.begin();
    CHECK(idx0 < idx1);
    CHECK(idx0 < idx2);
    CHECK(idx1 < idx3);
    CHECK(idx2 < idx3);
}

TEST_CASE("Poset topologically_sorted_list single element") {
    std::vector<std::pair<unsigned int, unsigned int>> empty_rel2;
    Poset p(1, empty_rel2);
    auto sorted = p.topologically_sorted_list();
    CHECK(sorted.size() == 1);
    CHECK(sorted[0] == 0);
}

// =============================================================================
// misc.h tests: read_objective
// =============================================================================

TEST_CASE("read_objective basic") {
    std::istringstream is("1.0 2.5 3.7 -0.5");
    std::vector<double> obj;
    read_objective(is, obj);
    CHECK(obj.size() == 4);
    CHECK(obj[0] == doctest::Approx(1.0));
    CHECK(obj[1] == doctest::Approx(2.5));
    CHECK(obj[2] == doctest::Approx(3.7));
    CHECK(obj[3] == doctest::Approx(-0.5));
}

TEST_CASE("read_objective empty stream") {
    std::istringstream is("");
    std::vector<double> obj;
    read_objective(is, obj);
    CHECK(obj.empty());
}

// =============================================================================
// misc.h tests: read_pointset
// =============================================================================

TEST_CASE("read_pointset basic numeric") {
    std::istringstream is("1 2 3\n4 5 6\n7 8 9\n");
    std::vector<std::vector<double>> pts;
    read_pointset(is, pts);
    CHECK(pts.size() == 3);
    CHECK(pts[0].size() == 3);
    CHECK(pts[0][0] == doctest::Approx(1));
    CHECK(pts[0][1] == doctest::Approx(2));
    CHECK(pts[0][2] == doctest::Approx(3));
    CHECK(pts[1][0] == doctest::Approx(4));
    CHECK(pts[2][2] == doctest::Approx(9));
}

TEST_CASE("read_pointset with fractions") {
    std::istringstream is("1/2 3/4\n2/3 5/6\n");
    std::vector<std::vector<double>> pts;
    read_pointset(is, pts);
    CHECK(pts.size() == 2);
    CHECK(pts[0][0] == doctest::Approx(0.5));
    CHECK(pts[0][1] == doctest::Approx(0.75));
}

TEST_CASE("read_pointset with negatives") {
    std::istringstream is("-1.5 2.0\n3.5 -4.0\n");
    std::vector<std::vector<double>> pts;
    read_pointset(is, pts);
    CHECK(pts.size() == 2);
    CHECK(pts[0][0] == doctest::Approx(-1.5));
    CHECK(pts[1][1] == doctest::Approx(-4.0));
}

TEST_CASE("read_pointset skip comment lines") {
    std::istringstream is("# this is a comment\n1 2\n# another comment\n3 4\n");
    std::vector<std::vector<double>> pts;
    read_pointset(is, pts);
    CHECK(pts.size() == 2);
}

TEST_CASE("read_pointset empty") {
    std::istringstream is("");
    std::vector<std::vector<double>> pts;
    read_pointset(is, pts);
    CHECK(pts.empty());
}

// =============================================================================
// misc.h tests: read_poset_from_file
// =============================================================================

TEST_CASE("read_poset_from_file basic") {
    std::istringstream is("3\n0 1\n0 2\n");
    Poset p = read_poset_from_file(is);
    CHECK(p.num_elem() == 3);
    CHECK(p.num_relations() == 2);
}

TEST_CASE("read_poset_from_file empty relations") {
    std::istringstream is("4\n");
    Poset p = read_poset_from_file(is);
    CHECK(p.num_elem() == 4);
    CHECK(p.num_relations() == 0);
}

// =============================================================================
// misc.h tests: read_poset_from_file_adj_matrix
// =============================================================================

TEST_CASE("read_poset_from_file_adj_matrix valid") {
    // 3x3 adjacency matrix: 0->1, 0->2, 1->2
    std::istringstream is("0 1 1\n0 0 1\n0 0 0\n");
    auto result = read_poset_from_file_adj_matrix(is);
    CHECK(result.first == true);
    CHECK(result.second.num_elem() == 3);
}

TEST_CASE("read_poset_from_file_adj_matrix invalid size") {
    // Matrix says n=3 but provides only 1x1
    std::istringstream is("1 1 1\n0\n");
    auto result = read_poset_from_file_adj_matrix(is);
    CHECK(result.first == false);
}

TEST_CASE("read_poset_from_file_adj_matrix empty") {
    std::istringstream is("0\n");
    auto result = read_poset_from_file_adj_matrix(is);
    CHECK(result.first == true);
    CHECK(result.second.num_elem() == 1);
}

// =============================================================================
// misc.h tests: read_inner_ball
// =============================================================================

#include "cartesian_geom/cartesian_kernel.h"

TEST_CASE("read_inner_ball basic") {
    std::istringstream is("1.0 2.0 0.5");
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    auto result = read_inner_ball<NT, Point>(is);
    CHECK(result.second == doctest::Approx(0.5));  // radius
}

// =============================================================================
// Poset is_in with int type
// =============================================================================

TEST_CASE("Poset is_in with int") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}};
    Poset p(2, rel);
    std::vector<int> pt1 = {1, 2};
    CHECK(p.is_in(pt1));
    std::vector<int> pt2 = {3, 2};
    CHECK(!p.is_in(pt2));
}

// =============================================================================
// Poset copy semantics
// =============================================================================

TEST_CASE("Poset copy constructor") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}};
    Poset p1(3, rel);
    Poset p2(p1);
    CHECK(p2.num_elem() == 3);
    CHECK(p2.num_relations() == 1);
}

TEST_CASE("Poset assignment") {
    std::vector<std::pair<unsigned int, unsigned int>> rel = {{0, 1}};
    Poset p1(3, rel);
    Poset p2;
    p2 = p1;
    CHECK(p2.num_elem() == 3);
    CHECK(p2.num_relations() == 1);
}
