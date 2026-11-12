// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis
// Copyright (c) 2021 Maios Papachristou

// Licensed under GNU LGPL.3, see LICENCE file


// Edited by HZ on 11.06.2020 - mute doctest.h
#include "doctest.h"
#include <fstream>
#include <iostream>



#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "misc/misc.h"
#include "random_walks/random_walks.hpp"

#include "volume/volume_sequence_of_balls.hpp"
#include "generators/known_polytope_generators.h"
#include "sampling/sampling.hpp"

#include "diagnostics/univariate_psrf.hpp"

#include "preprocess/inscribed_ellipsoid_rounding.hpp"


template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT get_samples(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    uniform_sampling<WalkType>(randPoints, P, rng, walkL, numpoints,
                               StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT get_samples_boundary(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    uniform_sampling_boundary<WalkType>(randPoints, P, rng, walkL, numpoints,
                                        StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

template
<
    typename MT,
    typename WalkType,
    typename Polytope
>
MT get_samples_gaussian(Polytope &P)
{
    typedef typename Polytope::PointType Point;
    typedef typename Polytope::NT NT;
    NT a = 1;

    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = P.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    gaussian_sampling<WalkType>(randPoints, P, rng, walkL, numpoints, a,
                               StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }

    return samples;
}

template <typename NT, typename WalkType = DikinWalk>
void call_test_dikin(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Dikin Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = JohnWalk>
void call_test_john(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing John Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = VaidyaWalk>
void call_test_vaidya(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Vaidya Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = BRDHRWalk>
void call_test_brdhr(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Boundary RDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_boundary<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = BCDHRWalk>
void call_test_bcdhr(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Boundary CDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_boundary<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = GaussianRDHRWalk>
void call_test_grdhr(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Gaussian RDHR for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_gaussian<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

// Every chord of a polytope far from the mode lies in the tail of the gaussian,
// where drawing from the 1-dimensional gaussian until it hits the chord used to
// practically never finish
template <typename NT, typename WalkType>
void call_test_gaussian_far_from_mode(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int d = 4, walkL = 10, numpoints = 10000, nburns = 0;
    NT a = 0.5; // the standard gaussian

    std::cout << "--- Testing gaussian sampling from the H-cube [9, 11]^4" << std::endl;
    Hpolytope P = generate_cube<Hpolytope>(d, false);
    P.shift(VT::Constant(d, NT(-10)));
    P.ComputeInnerBall();
    Point StartingPoint(d, std::vector<NT>(d, NT(10)));
    RNGType rng(d);
    std::list<Point> randPoints;

    gaussian_sampling<WalkType>(randPoints, P, rng, walkL, numpoints, a,
                                StartingPoint, nburns);

    MT samples(d, numpoints);
    unsigned int jj = 0, outside = 0;
    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit!=randPoints.end(); rpit++, jj++)
    {
        if (P.is_in(*rpit, NT(1e-10)) == 0) outside++;
        samples.col(jj) = (*rpit).getCoefficients();
    }
    CHECK(outside == 0);

    // the coordinates are independent, each one is the standard gaussian
    // truncated to [9, 11] with mean (phi(9) - phi(11)) / (Phi(11) - Phi(9))
    NT expected_mean = 9.108523101649192;
    VT mean = samples.rowwise().mean();
    std::cout << "max error of the mean = " << (mean.array() - expected_mean).abs().maxCoeff() << std::endl;
    CHECK((mean.array() - expected_mean).abs().maxCoeff() < 0.01);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;
    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = GaussianBallWalk>
void call_test_gbaw(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Gaussian Ball Walk for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_gaussian<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

template <typename NT, typename WalkType = GaussianHamiltonianMonteCarloExactWalk>
void call_test_ghmc(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing Gaussian exact HMC for H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    MT samples = get_samples_gaussian<MT, WalkType>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 2.2);
}

template <typename NT, typename WalkType = GaussianAcceleratedBilliardWalk>
void call_test_gabw(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0;
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    start = std::chrono::high_resolution_clock::now();

    std::cout << "--- Testing Gaussian Accelerated Billiard Walk for Skinny-H-cube10" << std::endl;
    P = generate_skinny_cube<Hpolytope>(d);


    Point p = P.ComputeInnerBall().first;
    typedef typename GaussianAcceleratedBilliardWalk::template Walk
            <
                    Hpolytope,
                    RNGType
            > walk;
    typedef MultivariateGaussianRandomPointGenerator <walk> RandomPointGenerator;
    PushBackWalkPolicy push_back_policy;

    EllipsoidParams<NT> params;
    params.john_params.maxiter = 500;
    params.john_params.tol = std::pow(10, -6.0);
    params.john_params.reg = std::pow(10, -4.0);

    std::tuple<MT, VT, NT> ellipsoid = compute_inscribed_ellipsoid<MT, EllipsoidType::MAX_ELLIPSOID>
    (P.get_mat(), P.get_vec(), p.getCoefficients(), params);
    const MT E = get<0>(ellipsoid);

    RandomPointGenerator::apply(P, p, E, numpoints, 1, randPoints,
                                push_back_policy, rng);

    stop = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> total_time = stop - start;
    std::cout << "Done in " << total_time.count() << '\n';

    MT samples(d, numpoints);
    unsigned int jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit != randPoints.end(); rpit++, jj++)
        samples.col(jj) = (*rpit).getCoefficients();

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
    RNGType Srng(d);

    typedef Eigen::SparseMatrix<NT> SparseMT;
    typedef HPolytope<Point, Eigen::SparseMatrix<NT, Eigen::RowMajor>> SparseHpoly;
    std::list<Point> Points;

    SparseHpoly SP;
    SP = generate_skinny_cube<SparseHpoly>(d);


    std::cout << "--- Testing Gaussian Accelerated Billiard Walk for Sparse Skinny-H-cube10" << std::endl;


    p = SP.ComputeInnerBall().first;
    typedef typename GaussianAcceleratedBilliardWalk::template Walk
            <
                    SparseHpoly,
                    RNGType,
                    SparseMT
            > sparsewalk;
    typedef MultivariateGaussianRandomPointGenerator <sparsewalk> SparseRandomPointGenerator;

    start = std::chrono::high_resolution_clock::now();

    ellipsoid = compute_inscribed_ellipsoid<MT, EllipsoidType::MAX_ELLIPSOID>
    ((SparseMT)SP.get_mat(), SP.get_vec(), p.getCoefficients(), params);

    const SparseMT SE = get<0>(ellipsoid).sparseView();

    SparseRandomPointGenerator::apply(SP, p, SE, numpoints, 1, Points,
                                push_back_policy, Srng);

    stop = std::chrono::high_resolution_clock::now();

    total_time = stop - start;
    std::cout << "Done in " << total_time.count() << '\n';

    jj = 0;
    MT Ssamples(d, numpoints);
    for (typename std::list<Point>::iterator rpit = Points.begin(); rpit != Points.end(); rpit++, jj++)
        Ssamples.col(jj) = (*rpit).getCoefficients();

    score = univariate_psrf<NT, VT>(Ssamples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
    NT delta = (samples - Ssamples).maxCoeff();
    std::cout << "\ndelta = " << delta << std::endl;
    CHECK(delta < 0.00001);

}

template <typename NT>
void call_test_sparse(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point, Eigen::SparseMatrix<NT, Eigen::RowMajor>> Hpolytope;
    typedef typename Hpolytope::MT MT;
    typedef typename Hpolytope::DenseMT DenseMT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    Hpolytope P;
    unsigned int d = 10;

    std::cout << "--- Testing ABW for sparse H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    DenseMT samples = get_samples<DenseMT, AcceleratedBilliardWalk>(P);

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);


    std::cout << "--- Testing CDHR for sparse H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    samples = get_samples<DenseMT, CDHRWalk>(P);

    score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);


    std::cout << "--- Testing RDHR for sparse H-cube10" << std::endl;
    P = generate_cube<Hpolytope>(d, false);
    P.ComputeInnerBall();

    samples = get_samples<DenseMT, RDHRWalk>(P);

    score = univariate_psrf<NT, VT>(samples);
    std::cout << "psrf = " << score.maxCoeff() << std::endl;

    CHECK(score.maxCoeff() < 1.1);
}

TEST_CASE("dikin") {
    call_test_dikin<double>();
}

TEST_CASE("john") {
    call_test_john<double>();
}

TEST_CASE("vaidya") {
    call_test_vaidya<double>();
}

TEST_CASE("brdhr") {
    call_test_brdhr<double>();
}

TEST_CASE("bcdhr") {
    call_test_bcdhr<double>();
}

TEST_CASE("grdhr") {
    call_test_grdhr<double>();
}

TEST_CASE("gaussian_far_from_mode") {
    call_test_gaussian_far_from_mode<double, GaussianCDHRWalk>();
    call_test_gaussian_far_from_mode<double, GaussianRDHRWalk>();
}

TEST_CASE("gbaw") {
    call_test_gbaw<double>();
}

TEST_CASE("ghmc") {
    call_test_ghmc<double>();
}

TEST_CASE("gabw") {
    call_test_gabw<double>();
}

TEST_CASE("sparse") {
    call_test_sparse<double>();
}
