#ifndef EIGCORRELATION
    #define EIGCORRELATION
#endif

#include <vector>
#include <chrono>
#include <iostream>
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"

#include "sampling/sample_correlation_matrices.hpp"

#include <fstream>

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
}

template <typename NT, typename WalkType, typename RNGType>
void naive_uniform_test(int n, int num_points, int walkL){
    typedef Cartesian<NT>                                               Kernel;
    typedef typename Kernel::Point                                      Point;
    typedef std::vector<Point>                                          PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 
    
    std::chrono::steady_clock::time_point start, end;
    double time;

    PointList randPoints;
    std::cout << "Direct implementation : " << n << " - time : ";
    
    start = std::chrono::steady_clock::now();

    direct_uniform_sampling<NT, WalkType, RNGType, Point>(n, num_points, walkL, randPoints, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << " (ms)" << std::endl;

    // write_to_file("sampling.txt", randPoints);
    // check_output<NT,VT,MT>(randPoints, num_points, n);
}

template <typename NT, typename WalkType, typename RNGType>
void new_test(unsigned int n, unsigned int const num_points, unsigned int walkL){

    std::cout << "Improved implementation : " << n << " - time : ";

    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef std::vector<Point> PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 

    std::chrono::steady_clock::time_point start, end;
    double time;

    PointList randPoints;

    start = std::chrono::steady_clock::now();
    
    uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);
    
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << " (ms)" << std::endl;

    // write_to_file("sampling_new.txt", randPoints);
    // check_output<NT,VT,MT>(randPoints, num_points, n);
}

template <typename NT, typename WalkType, typename RNGType>
void gaussian_test(unsigned int n, unsigned int const num_points, unsigned int walkL){

    std::cout << "Gaussian implementation : " << n << " - time : ";

    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef std::vector<Point> PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 
    
    std::chrono::steady_clock::time_point start, end;
    double time;

    PointList randPoints;
    
    start = std::chrono::steady_clock::now();
    
    gaussian_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 2);
    
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << " (ms)" << std::endl;

    // write_to_file("sampling_new.txt", randPoints);
    // check_output<NT,VT,MT>(randPoints, num_points, n);
}

int main(int argc, char const *argv[]) {
    srand((unsigned) time(NULL));
    // srand(19031999);
    typedef double NT;
    typedef Cartesian<NT>           Kernel;
    typedef typename Kernel::Point  Point;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT; 
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    
    // BilliardWalk, AcceleratedBilliardWalk, GaussianAcceleratedBilliardWalk
    // GaussianHamiltonianMonteCarloExactWalk
    unsigned int n, num_points = 100, walkL = 10;
    std::cout << "Input n = ";
    std::cin >> n;
    
    new_test<NT, BilliardWalk, RNGType>(n, num_points, walkL);

    naive_uniform_test<NT, BilliardWalk, RNGType>(n, num_points, walkL);

    // gaussian_test<NT, RNGType>(n, num_points, walkL);

    // naive_test<NT, GaussianHamiltonianMonteCarloExactWalk, RNGType>(n, num_points, walkL, nreflex);
    
    return 0;
}