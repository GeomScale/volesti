// #ifndefine EIGEN_USE_MKL_ALL
// #define EIGEN_USE_MKL_ALL

#include "doctest.h"
#include <fstream>
#include <iostream>
#include "misc.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"

#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sample_correlation_matrices.hpp"
#include "convex_bodies/correlation_matrices/corre_spectra2.hpp"

#include "diagnostics/univariate_psrf.hpp"

template <typename Point>
void write_to_file(std::string filename, std::vector<Point> const& randPoints) {
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf()); //save and redirect
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf); //reset to standard output again
}

template<typename NT, typename MT, typename VT>
MT rebuildMatrix(const VT &xvector, const unsigned int n){
    MT mat = MT::Identity(n,n);
    NT coeff;
    int i, j, ind = 0;
    for(i = 0; i < n ; ++i){
        for(j = i+1; j < n; ++j){
            coeff = xvector[ind];
            mat(i,j) = mat(j,i) = coeff;
            ++ind;
        }
    }
    return mat;
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

    if(num_points >= 100){
        MT samples(d, num_points);
        unsigned int jj = 0;

        for (typename PointList::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++){
            samples.col(jj) = (*rpit).getCoefficients();
        }

        VT score = univariate_psrf<NT, VT>(samples);
        std::cout << "psrf = " << score.maxCoeff() << std::endl;

        CHECK(score.maxCoeff() < 1.1);
    }
}

template<typename NT, typename VT, typename MT, typename PointList>
void check_output2(PointList &randPoints, int num_points, int n){
    int d = n*(n-1)/2, count = 0;
    int i, j, k;
    MT A;
    Eigen::LDLT<MT> A_ldlt;
    for(i = 0; i < num_points ; ++i){
        A = MT::Identity(n,n);
        for(j = 0; j < n ; ++j){
            for(k = j+1 ; k < n ; ++k){
                A(j,k) = A(k,j) = randPoints[i].mat(k,j);
            }
        }
        A_ldlt = Eigen::LDLT<MT>(A);
        if (A_ldlt.info() == Eigen::NumericalIssue || !A_ldlt.isPositive()){
            ++count;
        }
    }
    std::cout << "Fails " << count << " / " << num_points << " samples\n";
    CHECK(count == 0);

    // MT samples(d, num_points);
    // unsigned int jj = 0;

    // for (typename PointList::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++){
    //     samples.col(jj) = (*rpit).getCoefficients();
    // }

    // VT score = univariate_psrf<NT, VT>(samples);
    // std::cout << "psrf = " << score.maxCoeff() << std::endl;

    // CHECK(score.maxCoeff() < 1.1);
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

    CorreSpectra<Point> P(n); 
    Point startingPoint(n);

    CHECK(P.matrixSize() == n);
    CHECK(P.dimension() == d);
    CHECK(P.is_in(startingPoint) == -1);
    std::cout << "Diameter of P = " << P.InnerBall().second <<std::endl;

    CorreSpectra2<PMT> P2(n); 

    CHECK(P2.matrixSize() == n);
    CHECK(P2.dimension() == d);

    PMT startingPoint2(n);
    // compute_diameter<CorreSpectra2<PMT>>::compute<NT>(P2);
    PMT A = GetDirection<PMT>::apply(P2.dimension(), rng);
}

template
<
    typename NT,
    typename WalkType
>
void test_old_uniform(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>           RNGType;

    std::cout << "Test old sampling : "<< num_points <<" uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;


    start = std::chrono::steady_clock::now();

    direct_uniform_sampling<NT, WalkType, RNGType, Point>(n, num_points, walkL, randPoints, 0);

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
void test_new_uniform2(const unsigned int n, const unsigned int num_points = 1000){
    typedef CorreMatrix<NT>                                     Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;
    std::cout << "Test new sampling 2 : "<< num_points << " uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();
    
    uniform_correlation_sampling2<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);
    
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    // check_output2<NT, VT, MT>(randPoints, num_points, n);
}

template 
<
    typename NT, 
    typename WalkType
>
void test_new_gaussian(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

    std::cout << "Test new sampling : "<< num_points << " gaussian correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    // gaussian_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, NT(0.5));

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output<NT, VT, MT>(randPoints, num_points, n);
}

TEST_CASE("corre_spectra") {
    test_corre_spectra_classes<double>(10);
}

///////////////////////////////////////////////////////////////////
//      Old implementation
//

// TEST_CASE("old_ball_uniform") {
//     test_old_uniform<double, BallWalk>(10,100);
// }

TEST_CASE("old_billiard_uniform") {
    test_old_uniform<double, BilliardWalk>(10,100);
}

TEST_CASE("old_accelerated_billiard_uniform") {
    test_old_uniform<double, AcceleratedBilliardWalk>(10,100);
}

///////////////////////////////////////////////////////////////////
//      New implementation : CorreSpectra Vector PointType
//

TEST_CASE("new_ball_uniform") {
    test_new_uniform<double, BallWalk>(10,1000000);
}

TEST_CASE("new_billiard_uniform") {
    test_new_uniform<double, BilliardWalk>(10,100);
}

TEST_CASE("new_accelerated_billiard_uniform") {
    test_new_uniform<double, AcceleratedBilliardWalk>(3,1);
}

///////////////////////////////////////////////////////////////////
//      New implementation : CorreSpectra Matrix PointType
//

TEST_CASE("new_billiard_uniform_2") {
    test_new_uniform2<double, BilliardWalk>(3,1);
}

TEST_CASE("new_accelerated_billiard_uniform_2") {
    test_new_uniform2<double, AcceleratedBilliardWalk>(10,100);
}

// TEST_CASE("new_ReHMC_gaussian") {
//     call_test_new_ReHMC_gaussian<double>(3);
// }

// #endif