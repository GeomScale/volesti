// Licensed under GNU LGPL.3, see LICENCE file
// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2021 Vissarion Fisikopoulos
// Copyright (c) 2018-2021 Apostolos Chalkis

// Contributed by Huu Phuoc Le as part of Google Summer of Code 2022 program

// Licensed under GNU LGPL.3, see LICENCE file

#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sample_correlation_matrices.hpp"

typedef double                                              NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>   MT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1>                VT;

typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3>   RNGType;

typedef Cartesian<NT>                                       Kernel;
typedef typename Kernel::Point                              Point;
typedef CorreMatrix<NT>                                     PointMT;

// Direct call to sampling algorithms for Spectrahedron class for comparison

std::vector<MT> LMIGen(int n){
    int i, j, l, k = n*(n-1)/2+1;
    std::vector<MT> list_Mat;
    MT A;
    list_Mat.push_back(-MT::Identity(n, n));
    for(i = 0; i < n; ++i){
        for(j = 0; j < i; ++j){
            A = MT::Zero(n, n);
            A(i,j) = -1;
            A(j,i) = -1;
            list_Mat.push_back(A);
        }
    }
    return list_Mat;
}

Spectrahedron<Point> prepare_input(int n){
    int d = n*(n-1)/2;
    Point p(d);
    std::vector<MT> lmi_mat = LMIGen(n);
    LMI<NT, MT, VT> lmi(lmi_mat);
    Spectrahedron<Point> spectra(lmi);
    spectra.set_interior_point(p);
    spectra._inner_ball.second = 1/std::sqrt(d);
    return spectra;
}

template<typename WalkType>
void old_uniform_sampling(const unsigned int n, const unsigned int num_points){

    std::cout << "Old sampling : "<< num_points <<" uniform correlation matrices of size " << n << std::endl;
    std::chrono::steady_clock::time_point start, end;
    double time;
    std::vector<Point> randPoints;
    unsigned int walkL = 1;

    Spectrahedron<Point> spectra = prepare_input(n);
    int d = spectra.dimension();
    Point p(d);
    RNGType rng(d);

    start = std::chrono::steady_clock::now();

    uniform_sampling<WalkType>(randPoints, spectra, rng, walkL, num_points, p, 0);

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;
}

template<typename PointType>
void write_to_file(std::string filename, std::vector<PointType> const& randPoints){
    std::ofstream out(filename);
    auto coutbuf = std::cout.rdbuf(out.rdbuf());
    for(int i=0; i<randPoints.size(); ++i)
        randPoints[i].print();
    std::cout.rdbuf(coutbuf);
}

template<typename WalkType>
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

template<typename WalkType>
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

int main(int argc, char const *argv[]){

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

    unsigned int n = 3, num_points = 5000;

    old_uniform_sampling<BilliardWalk>(n, num_points);

    correlation_matrix_uniform_sampling<BallWalk>(n, num_points, "BallWalk");

    correlation_matrix_uniform_sampling<RDHRWalk>(n, num_points, "RDHRWalk");

    correlation_matrix_uniform_sampling<BilliardWalk>(n, num_points, "BilliardWalk");

    correlation_matrix_uniform_sampling<AcceleratedBilliardWalk>(n, num_points, "AcceleratedBilliardWalk");

    correlation_matrix_uniform_sampling_MT<BallWalk>(n, num_points, "BallWalk");

    correlation_matrix_uniform_sampling_MT<RDHRWalk>(n, num_points, "RDHRWalk");

    correlation_matrix_uniform_sampling_MT<BilliardWalk>(n, num_points, "BilliardWalk");

    correlation_matrix_uniform_sampling_MT<AcceleratedBilliardWalk>(n, num_points, "AcceleratedBilliardWalk");

    return 0;
}