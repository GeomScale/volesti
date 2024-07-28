// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <iostream>

#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "cartesian_geom/cartesian_kernel.h"
#include "diagnostics/univariate_psrf.hpp"
#include "misc/misc.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sample_correlation_matrices.hpp"

template<typename NT, typename MT, typename VT>
MT rebuildMatrix(const VT &xvector, const unsigned int n){
    MT mat = MT::Identity(n,n);
    NT coeff;
    int i, j, ind = 0;
    for(i = 0; i < n ; ++i){
        for(j = 0; j < i; ++j){
            coeff = xvector[ind];
            mat(i,j) = mat(j,i) = coeff;
            ++ind;
        }
    }
    return mat;
}

template<typename NT, typename MT>
Eigen::Matrix<NT, Eigen::Dynamic, 1> getCoefficientsFromMatrix(const MT& mat) {
    int n = mat.rows();
    int d = n * (n - 1) / 2;
    Eigen::Matrix<NT, Eigen::Dynamic, 1> coeffs(d);
    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            coeffs(k) = mat(i, j);
            ++k;
    	}
    }
    return coeffs;
}

template<typename NT, typename VT, typename MT, typename PointList>
void check_output(PointList &randPoints, int num_points, int n){
    int d = n*(n-1)/2, count = 0;
    MT A;
    Eigen::LDLT<MT> A_ldlt;
    for(int i = 0; i < num_points ; ++i){
        A = rebuildMatrix<NT, MT>(randPoints[i].getCoefficients(), n);
        A_ldlt = Eigen::LDLT<MT>(A);
        if (A_ldlt.info() == Eigen::NumericalIssue || !A_ldlt.isPositive()){
            ++count;
        }
    }
    std::cout << "Fails " << count << " / " << num_points << " samples\n";
    CHECK(count == 0);

    MT samples(d, num_points);
    unsigned int jj = 0;

    for (typename PointList::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++){
        samples.col(jj) = (*rpit).getCoefficients();
    }

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template<typename NT, typename VT, typename MT>
void check_output_MT(std::list<MT> &randCorMatrices, int num_points, int n){
    int d = n*(n-1)/2, count = 0;
    MT A;
    Eigen::LDLT<MT> mat_ldlt;
    for(auto& mat : randCorMatrices){
        mat_ldlt = Eigen::LDLT<MT>(mat);
    	if(mat_ldlt.info() == Eigen::NumericalIssue || !mat_ldlt.isPositive()){
            ++count;
    	}
    }
    std::cout << "Fails " << count << " / " << num_points << " samples\n";
    CHECK(count == 0);

    MT samples(d, num_points);
    unsigned int jj = 0;
    for(const auto& mat : randCorMatrices){
        samples.col(jj) = getCoefficientsFromMatrix<NT, MT>(mat);
    	jj++;
    }

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT>
void test_corre_spectra_classes(unsigned int const n){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>     MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                  VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;
    typedef CorreMatrix<NT>                                     PMT;

    const unsigned int d = n*(n-1)/2;
    RNGType rng(d);

    CorrelationSpectrahedron<Point> P(n);
    Point startingPoint(d);

    CHECK(P.matrixSize() == n);
    CHECK(P.dimension() == d);
    CHECK(P.is_in(startingPoint) == -1);
    std::cout << "Diameter of P = " << P.InnerBall().second <<std::endl;

    CorrelationSpectrahedron_MT<PMT> P2(n);

    CHECK(P2.matrixSize() == n);
    CHECK(P2.dimension() == d);

    PMT startingPoint2(d);
    CHECK(P2.is_in(startingPoint2) == -1);
    PMT A = GetDirection<PMT>::apply(P2.dimension(), rng);
}

template
<
    typename NT,
    typename WalkType
>
void test_new_uniform(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

    std::cout << "Test new sampling : "<< num_points << " uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output<NT, VT, MT>(randPoints, num_points, n);
}

template
<
    typename NT,
    typename WalkType
>
void test_new_uniform_MT(const unsigned int n, const unsigned int num_points = 1000){
    typedef CorreMatrix<NT>                                     Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;
    std::cout << "Test new sampling 2 : "<< num_points << " uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::list<MT> randCorMatrices;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    uniform_correlation_sampling_MT<WalkType, Point, RNGType>(n, randCorMatrices, walkL, num_points, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output_MT<NT, VT, MT>(randCorMatrices, num_points, n);
}

int n = 3;
int num_points_BallWalk = 10000;
int num_points_RDHRWalk = 10000;
int num_points_BilliardWalk = 1000;

///////////////////////////////////////////////////////////////////
//      Test new classes
//

TEST_CASE("corre_spectra") {
    test_corre_spectra_classes<double>(n);
}

///////////////////////////////////////////////////////////////////
//      New implementation : CorrelationSpectrahedron Vector PointType
//

TEST_CASE("new_ball_uniform") {
    std::cout << "Ball Walk :: ";
    test_new_uniform<double, BallWalk>(n,num_points_BallWalk);
}

TEST_CASE("new_billiard_uniform") {
    std::cout << "Billiard Walk :: ";
    test_new_uniform<double, BilliardWalk>(n, num_points_BilliardWalk);
}

TEST_CASE("new_accelerated_billiard_uniform") {
    std::cout << "Accelerated Billiard Walk :: ";
    test_new_uniform<double, AcceleratedBilliardWalk>(n, num_points_BilliardWalk);
}

///////////////////////////////////////////////////////////////////
//      New implementation : CorrelationSpectrahedron Matrix PointType
//

TEST_CASE("new_ball_uniform_MT") {
    std::cout << "Ball Walk MT :: ";
    test_new_uniform_MT<double, BallWalk>(n,num_points_BallWalk);
}

TEST_CASE("new_rdhr_uniform_MT") {
    std::cout << "RDHR Walk MT :: ";
    test_new_uniform_MT<double, RDHRWalk>(n,num_points_RDHRWalk);
}

TEST_CASE("new_billiard_uniform_MT") {
    std::cout << "Billiard Walk MT :: ";
    test_new_uniform_MT<double, BilliardWalk>(n, num_points_BilliardWalk);
}

TEST_CASE("new_accelerated_billiard_uniform_MT") {
    std::cout << "Accelerated Billiard Walk MT :: ";
    test_new_uniform_MT<double, AcceleratedBilliardWalk>(n, num_points_BilliardWalk);
}
