#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sample_correlation_matrices.hpp"

template <typename Point>
void write_to_file(std::string filename, std::vector<Point> const& randPoints) {
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf());
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf);
}

template
<
    typename NT,
    typename WalkType
>
void correlation_matrix_uniform_sampling(const unsigned int n, const unsigned int num_points = 1000){
    typedef Cartesian<NT>                                       Kernel;
    typedef typename Kernel::Point                              Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

    std::cout << "Sample uniformly "<< num_points << " correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<Point>("correlation_matrices.txt", randPoints);
}

template
<
    typename NT,
    typename WalkType
>
void correlation_matrix_uniform_sampling_MT(const unsigned int n, const unsigned int num_points = 1000){
    typedef CorreMatrix<NT>                                     Point;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

    std::cout << "Sample uniformly "<< num_points << " correlation matrices of size " << n << " with matrix PointType" << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();
    
    uniform_correlation_sampling_MT<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);
    
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<Point>("correlation_matrices_MT.txt", randPoints);
}

int main(int argc, char const *argv[]) {
    srand((unsigned) time(NULL));
    typedef double NT;
    unsigned int n = 3, num_points = 5000;
    
    // correlation_matrix_uniform_sampling<NT, BilliardWalk>(n,num_points);

    correlation_matrix_uniform_sampling<NT, AcceleratedBilliardWalk>(n,num_points);
    
    // correlation_matrix_uniform_sampling_MT<NT, BilliardWalk>(n,num_points);

    // correlation_matrix_uniform_sampling_MT<NT, AcceleratedBilliardWalk>(n,num_points);
    
    return 0;
}