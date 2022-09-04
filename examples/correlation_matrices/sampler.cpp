#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sample_correlation_matrices.hpp"

typedef double                                              NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

typedef Cartesian<NT>                                       Kernel;
typedef typename Kernel::Point                              Point;

typedef CorreMatrix<NT>                                     PointMT;

template <typename PointType>
void write_to_file(std::string filename, std::vector<PointType> const& randPoints) {
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf());
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf);
}

template <typename WalkType>
void correlation_matrix_uniform_sampling(const unsigned int n, const unsigned int num_points, std::string walkname){

    std::cout << walkname << " samples uniformly "<< num_points << " correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<Point>(walkname + "_matrices.txt", randPoints);
}

template <typename WalkType>
void correlation_matrix_uniform_sampling_MT(const unsigned int n, const unsigned int num_points, std::string walkname){
    
    std::cout << walkname << " samples uniformly "<< num_points << " correlation matrices of size " << n << " with matrix PointType" << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<PointMT> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();
    
    uniform_correlation_sampling_MT<WalkType, PointMT, RNGType>(n, randPoints, walkL, num_points, 0);
    
    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<PointMT>(walkname + "_matrices_MT.txt", randPoints);
}

template <typename WalkType>
void correlation_matrix_gaussian_sampling(const unsigned int n, const unsigned int num_points, NT a, std::string walkname){

    std::cout << walkname << " samples "<< num_points << " correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    gaussian_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, a, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<Point>(walkname + "_matrices.txt", randPoints);
}

template <typename WalkType>
void correlation_matrix_gaussian_sampling_MT(const unsigned int n, const unsigned int num_points, NT a, std::string walkname){

    std::cout << walkname << " samples "<< num_points << " correlation matrices of size " << n << " with Matrix PointType" << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<PointMT> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    gaussian_correlation_sampling_MT<WalkType, PointMT, RNGType>(n, randPoints, walkL, num_points, a, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<PointMT>(walkname + "_matrices_MT.txt", randPoints);
}

template <typename WalkType>
void correlation_matrix_exponential_sampling(const unsigned int n, const unsigned int num_points, VT c, NT T, std::string walkname){

    std::cout << walkname << " samples "<< num_points << " correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    exponential_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, c, T, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<Point>(walkname + "_matrices.txt", randPoints);
}

template <typename WalkType>
void correlation_matrix_exponential_sampling_MT(const unsigned int n, const unsigned int num_points, VT c, NT T, std::string walkname){

    std::cout << walkname << " samples "<< num_points << " correlation matrices of size " << n << " with Matrix PointType" << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<PointMT> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    exponential_correlation_sampling_MT<WalkType, PointMT, RNGType>(n, randPoints, walkL, num_points, c, T, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    write_to_file<PointMT>(walkname + "_matrices_MT.txt", randPoints);
}

int main(int argc, char const *argv[]) {

// To enable Intel MKL, change option USE_MKL to ON in CMakeLists.txt 

#ifdef EIGEN_USE_MKL_ALL

    MKLVersion Version;
 
    mkl_get_version(&Version);
 
    printf("Using Intel MKL %d.%d.%d\n",Version.MajorVersion,Version.MinorVersion,Version.UpdateVersion);
    printf("Platform:                %s\n",Version.Platform);
    printf("Processor optimization:  %s\n",Version.Processor);
    printf("================================================================\n");
    printf("\n");
#endif

    unsigned int n = 3, d = n*(n-1)/2, num_points = 10000;
    
    ///////////////////////////////////////////////////////////
    //          Uniform sampling
    
    correlation_matrix_uniform_sampling<BallWalk>(n, num_points, "BallWalk");

    correlation_matrix_uniform_sampling<RDHRWalk>(n, num_points, "RDHRWalk");
    
    correlation_matrix_uniform_sampling<BilliardWalk>(n, num_points, "BilliardWalk");

    correlation_matrix_uniform_sampling<AcceleratedBilliardWalk>(n, num_points, "AcceleratedBilliardWalk");

    correlation_matrix_uniform_sampling_MT<BallWalk>(n, num_points, "BallWalk");

    correlation_matrix_uniform_sampling_MT<RDHRWalk>(n, num_points, "RDHRWalk");
    
    correlation_matrix_uniform_sampling_MT<BilliardWalk>(n, num_points, "BilliardWalk");

    correlation_matrix_uniform_sampling_MT<AcceleratedBilliardWalk>(n, num_points, "AcceleratedBilliardWalk");
    
    ///////////////////////////////////////////////////////////
    //          Gaussian sampling
    
    NT a = NT(2);

    correlation_matrix_gaussian_sampling<GaussianReHMCWalk>(n, num_points, a, "GuassianReHMC");
    correlation_matrix_gaussian_sampling_MT<GaussianReHMCWalk>(n, num_points, a, "GuassianReHMC");

    ///////////////////////////////////////////////////////////
    //          Exponential sampling

    VT c(d);
    for(int i = 0; i < d; ++i){
        c(i) = 1;
    }
    NT T = NT(1);

    correlation_matrix_exponential_sampling<ExponentialReHMCWalk>(n, num_points, c, T, "ExponentialReHMC");
    correlation_matrix_exponential_sampling_MT<ExponentialReHMCWalk>(n, num_points, c, T, "ExponentialReHMC");
    
    return 0;
}