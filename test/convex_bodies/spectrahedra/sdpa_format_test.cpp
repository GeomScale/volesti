// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2025 Vissarion Fisikopoulos
// Copyright (c) 2020-2025 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"

#include <fstream>
#include <cstdio>
#include <unistd.h>
#include <list>
#include <string>
#include <iostream>

#include "SDPAFormatManager.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "cartesian_geom/cartesian_kernel.h"

template <typename NT>
void compare_matrix_vectors(
    std::vector<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>> const& a,
    std::vector<Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>> const& b) {
    CHECK(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].rows() == b[i].rows());
        CHECK(a[i].cols() == b[i].cols());
        for (int r = 0; r < a[i].rows(); ++r) {
            for (int c = 0; c < a[i].cols(); ++c) {
                CHECK(a[i](r, c) == doctest::Approx(b[i](r, c)));
            }
        }
    }
}

TEST_CASE("sdpa_raw_roundtrip") {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    SdpaFormatManager<NT> manager;

    int const dim = 2;
    int const m = 2;

    // Build matrices:
    //   A0 = I
    //   A1 = [[0, 1], [1, 0]]
    //   A2 = [[2, 0], [0, 3]]
    std::vector<MT> original_matrices(dim + 1, MT::Zero(m, m));
    original_matrices[0] = MT::Identity(m, m);
    original_matrices[1] << 0, 1,
                            1, 0;
    original_matrices[2] << 2, 0,
                            0, 3;

    VT original_obj(dim);
    original_obj << 5.0, 6.0;

    // Temporary file
    char tmpname[] = "/tmp/sdpa_test_raw_XXXXXX";
    int fd = mkstemp(tmpname);
    REQUIRE(fd != -1);
    close(fd);

    // Write
    {
        std::ofstream ofs(tmpname);
        REQUIRE(ofs.is_open());
        manager.writeSDPAFormatFile(ofs, original_matrices, original_obj);
        ofs.close();
    }

    // Read back
    std::vector<MT> read_matrices;
    VT read_obj;
    {
        std::ifstream ifs(tmpname);
        REQUIRE(ifs.is_open());
        manager.loadSDPAFormatFile(ifs, read_matrices, read_obj);
        ifs.close();
    }

    // Compare
    compare_matrix_vectors(original_matrices, read_matrices);
    CHECK(original_obj.size() == read_obj.size());
    for (int i = 0; i < original_obj.size(); ++i) {
        CHECK(original_obj(i) == doctest::Approx(read_obj(i)));
    }

    std::remove(tmpname);
}

TEST_CASE("sdpa_spectrahedron_roundtrip") {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Cartesian<NT>::Point Point;

    SdpaFormatManager<NT> manager;

    int const dim = 2;
    int const m = 2;

    // Build matrices
    std::vector<MT> matrices(dim + 1, MT::Zero(m, m));
    matrices[0] = MT::Identity(m, m);
    matrices[1] << 0, 1,
                   1, 0;
    matrices[2] << 2, 0,
                   0, 3;

    // Build LMI and Spectrahedron
    LMI<NT, MT, VT> lmi(matrices);
    Spectrahedron<Point> spectrahedron(lmi);

    // Build objective function Point
    VT obj_coeffs(dim);
    obj_coeffs << 5.0, 6.0;
    Point original_obj(obj_coeffs);

    // Temporary file
    char tmpname[] = "/tmp/sdpa_test_spectra_XXXXXX";
    int fd = mkstemp(tmpname);
    REQUIRE(fd != -1);
    close(fd);

    // Write
    {
        std::ofstream ofs(tmpname);
        REQUIRE(ofs.is_open());
        manager.writeSDPAFormatFile(ofs, spectrahedron, original_obj);
        ofs.close();
    }

    // Read back
    Spectrahedron<Point> read_spectrahedron;
    Point read_obj;
    {
        std::ifstream ifs(tmpname);
        REQUIRE(ifs.is_open());
        manager.loadSDPAFormatFile(ifs, read_spectrahedron, read_obj);
        ifs.close();
    }

    // Compare matrices via LMI
    std::vector<MT> lmi_orig = spectrahedron.getLMI().getMatrices();
    std::vector<MT> lmi_read = read_spectrahedron.getLMI().getMatrices();
    compare_matrix_vectors(lmi_orig, lmi_read);

    // Compare objective function coefficients
    VT orig_coeffs = original_obj.getCoefficients();
    VT read_coeffs = read_obj.getCoefficients();
    CHECK(orig_coeffs.size() == read_coeffs.size());
    for (int i = 0; i < orig_coeffs.size(); ++i) {
        CHECK(orig_coeffs(i) == doctest::Approx(read_coeffs(i)));
    }

    std::remove(tmpname);
}

TEST_CASE("sdpa_comment_lines") {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    SdpaFormatManager<NT> manager;

    // Manually construct an SDPA file with comment lines
    char tmpname[] = "/tmp/sdpa_test_comments_XXXXXX";
    int fd = mkstemp(tmpname);
    REQUIRE(fd != -1);
    close(fd);

    {
        std::ofstream ofs(tmpname);
        REQUIRE(ofs.is_open());
        ofs << "\" This is a comment line\n";
        ofs << "* Another comment\n";
        ofs << "2\n";
        ofs << "1\n";
        ofs << "2\n";
        ofs << "5.0 6.0\n";
        ofs << "1 0\n";
        ofs << "0 1\n";
        ofs << "0 -1\n";
        ofs << "-1 0\n";
        ofs << "-2 0\n";
        ofs << "0 -3\n";
        ofs.close();
    }

    std::vector<MT> read_matrices;
    VT read_obj;
    {
        std::ifstream ifs(tmpname);
        REQUIRE(ifs.is_open());
        manager.loadSDPAFormatFile(ifs, read_matrices, read_obj);
        ifs.close();
    }

    // Verify structure
    CHECK(read_matrices.size() == 3);  // A0, A1, A2
    CHECK(read_obj.size() == 2);
    CHECK(read_obj(0) == doctest::Approx(5.0));
    CHECK(read_obj(1) == doctest::Approx(6.0));

    // A0 should be identity
    CHECK(read_matrices[0](0, 0) == doctest::Approx(1.0));
    CHECK(read_matrices[0](0, 1) == doctest::Approx(0.0));
    CHECK(read_matrices[0](1, 0) == doctest::Approx(0.0));
    CHECK(read_matrices[0](1, 1) == doctest::Approx(1.0));

    // A1 should be [[0, 1], [1, 0]] — the file stores -A1, so read negates back
    CHECK(read_matrices[1](0, 0) == doctest::Approx(0.0));
    CHECK(read_matrices[1](0, 1) == doctest::Approx(1.0));
    CHECK(read_matrices[1](1, 0) == doctest::Approx(1.0));
    CHECK(read_matrices[1](1, 1) == doctest::Approx(0.0));

    // A2 should be [[2, 0], [0, 3]]
    CHECK(read_matrices[2](0, 0) == doctest::Approx(2.0));
    CHECK(read_matrices[2](0, 1) == doctest::Approx(0.0));
    CHECK(read_matrices[2](1, 0) == doctest::Approx(0.0));
    CHECK(read_matrices[2](1, 1) == doctest::Approx(3.0));

    std::remove(tmpname);
}

TEST_CASE("sdpa_diagonal_blocks") {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    SdpaFormatManager<NT> manager;

    // File layout:
    //   1 variable, block structure {2, -1} (2x2 dense + 1x1 diagonal)
    //   total matrix dimension = 3
    //   Matrices are read with blockStructure {2, -1}, meaning:
    //     - First block: 2x2 dense (4 values, row-major)
    //     - Second block: 1x1 diagonal (1 value)
    //
    // A0 (stored as-is):  [[1, 0, 0],
    //                      [0, 1, 0],
    //                      [0, 0, 3]]
    // A1 file lines store -A1 so that after negation we get:
    //                     [[4, 1, 0],
    //                      [1, 5, 0],
    //                      [0, 0, 6]]
    char tmpname[] = "/tmp/sdpa_test_diag_XXXXXX";
    int fd = mkstemp(tmpname);
    REQUIRE(fd != -1);
    close(fd);

    {
        std::ofstream ofs(tmpname);
        REQUIRE(ofs.is_open());
        ofs << "1\n";
        ofs << "2\n";
        ofs << "2 -1\n";
        ofs << "42.0\n";
        // A0 dense block rows
        ofs << "1 0\n";
        ofs << "0 1\n";
        // A0 diagonal block
        ofs << "3\n";
        // -A1 dense block rows (these will be negated on read to produce A1)
        ofs << "-4 -1\n";
        ofs << "-1 -5\n";
        // -A1 diagonal block
        ofs << "-6\n";
        ofs.close();
    }

    std::vector<MT> read_matrices;
    VT read_obj;
    {
        std::ifstream ifs(tmpname);
        REQUIRE(ifs.is_open());
        manager.loadSDPAFormatFile(ifs, read_matrices, read_obj);
        ifs.close();
    }

    // Verify
    CHECK(read_matrices.size() == 2);  // A0 and A1
    CHECK(read_obj.size() == 1);
    CHECK(read_obj(0) == doctest::Approx(42.0));

    // A0 — stored as-is
    CHECK(read_matrices[0](0, 0) == doctest::Approx(1.0));
    CHECK(read_matrices[0](0, 1) == doctest::Approx(0.0));
    CHECK(read_matrices[0](0, 2) == doctest::Approx(0.0));
    CHECK(read_matrices[0](1, 0) == doctest::Approx(0.0));
    CHECK(read_matrices[0](1, 1) == doctest::Approx(1.0));
    CHECK(read_matrices[0](1, 2) == doctest::Approx(0.0));
    CHECK(read_matrices[0](2, 0) == doctest::Approx(0.0));
    CHECK(read_matrices[0](2, 1) == doctest::Approx(0.0));
    CHECK(read_matrices[0](2, 2) == doctest::Approx(3.0));

    // A1 — read from file and negated
    CHECK(read_matrices[1](0, 0) == doctest::Approx(4.0));
    CHECK(read_matrices[1](0, 1) == doctest::Approx(1.0));
    CHECK(read_matrices[1](0, 2) == doctest::Approx(0.0));
    CHECK(read_matrices[1](1, 0) == doctest::Approx(1.0));
    CHECK(read_matrices[1](1, 1) == doctest::Approx(5.0));
    CHECK(read_matrices[1](1, 2) == doctest::Approx(0.0));
    CHECK(read_matrices[1](2, 0) == doctest::Approx(0.0));
    CHECK(read_matrices[1](2, 1) == doctest::Approx(0.0));
    CHECK(read_matrices[1](2, 2) == doctest::Approx(6.0));

    std::remove(tmpname);
}

TEST_CASE("sdpa_error_truncated_file") {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    SdpaFormatManager<NT> manager;

    // File with only the header, no matrix data
    char tmpname[] = "/tmp/sdpa_test_trunc_XXXXXX";
    int fd = mkstemp(tmpname);
    REQUIRE(fd != -1);
    close(fd);

    {
        std::ofstream ofs(tmpname);
        REQUIRE(ofs.is_open());
        ofs << "2\n1\n2\n1.0 2.0\n";
        ofs.close();
    }

    std::vector<MT> matrices;
    VT obj;
    {
        std::ifstream ifs(tmpname);
        REQUIRE(ifs.is_open());
        CHECK_THROWS_AS(manager.loadSDPAFormatFile(ifs, matrices, obj), int);
        ifs.close();
    }

    std::remove(tmpname);
}
