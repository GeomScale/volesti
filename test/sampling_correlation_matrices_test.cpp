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

    CHECK(P.matrixSize() == n);
    CHECK(P.dimension() == d);

    CorreSpectra2<PMT> P2(n); 

    CHECK(P2.matrixSize() == n);
    CHECK(P2.dimension() == d);

    PMT startingPoint(n);
    // compute_diameter<CorreSpectra2<PMT>>::compute<NT>(P2);
    PMT A = GetDirection<PMT>::apply(P2.dimension(), rng);
}

template
<
    typename NT,
    typename WalkType,
    typename PointList
>
void test_old_uniform_correlation_matrices(unsigned int n, const unsigned int num_points, PointList &randPoints)
{
    typedef Cartesian<NT>                                               Kernel;
    typedef typename Kernel::Point                                      Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>           RNGType;

    unsigned int walkL = 1;

    direct_uniform_sampling<NT, WalkType, RNGType, Point>(n, num_points, walkL, randPoints, 0);
}

template
<
    typename NT,
    typename WalkType,
    typename PointList
>
void test_new_uniform_correlation_matrices(unsigned int n, const unsigned int num_points, PointList &randPoints)
{
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

    unsigned int walkL = 1;

    uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);
}

template
<
    typename NT,
    typename WalkType,
    typename PointList
>
void test_new_uniform_correlation_matrices2(unsigned int n, const unsigned int num_points, PointList &randPoints)
{
    typedef CorreMatrix<NT>                                     Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

    unsigned int walkL = 1;

    uniform_correlation_sampling2<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);
}

template 
<
    typename NT, 
    typename WalkType,
    typename PointList
>
void test_new_gaussian_correlation_matrices(unsigned int n, const unsigned int num_points, PointList &randPoints)
{
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;
    
    unsigned int walkL = 1;

    gaussian_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, NT(0.5));
}

template<typename NT>
void call_test_old_billiard(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    
    std::cout << "Test old sampling : 1000 uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;

    start = std::chrono::steady_clock::now();

    test_old_uniform_correlation_matrices<double, BilliardWalk>(n, num_points, randPoints);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output<NT, VT, MT>(randPoints, num_points, n);
}

template<typename NT>
void call_test_new_billiard(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

    std::cout << "Test new sampling : 1000 uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;

    start = std::chrono::steady_clock::now();

    test_new_uniform_correlation_matrices<double, BilliardWalk>(n, num_points, randPoints);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output<NT, VT, MT>(randPoints, num_points, n);
}

template<typename NT>
void call_test_new_billiard2(const unsigned int n, const unsigned int num_points = 1000){
    typedef CorreMatrix<NT>                                     Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

    std::cout << "Test new sampling 2 : 1000 uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;

    start = std::chrono::steady_clock::now();

    test_new_uniform_correlation_matrices2<double, BilliardWalk>(n, num_points, randPoints);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output2<NT, VT, MT>(randPoints, num_points, n);
}

template<typename NT>
void call_test_new_ReHMC_gaussian(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

    std::cout << "Test new sampling : 1000 gaussian correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;

    start = std::chrono::steady_clock::now();

    // test_new_gaussian_correlation_matrices<double, GaussianReHMCCorrelationWalk>(n, randPoints);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    check_output<NT, VT, MT>(randPoints, num_points, n);
}

TEST_CASE("corre_spectra") {
    test_corre_spectra_classes<double>(10);
}

TEST_CASE("old_billiard_uniform") {
    call_test_old_billiard<double>(10,100);
}

TEST_CASE("new_billiard_uniform") {
    call_test_new_billiard<double>(10,100);
}

TEST_CASE("new_billiard_uniform_2") {
    call_test_new_billiard2<double>(10,100);
}

// TEST_CASE("new_ReHMC_gaussian") {
//     call_test_new_ReHMC_gaussian<double>(3);
// }
