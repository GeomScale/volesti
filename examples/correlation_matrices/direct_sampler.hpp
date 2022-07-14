// Direct implementation: call to sampling algorithms for spectrahedra in volesti
#include "sampling/sampling.hpp"
#include "auxiliaries.h"

template <typename MT>
std::vector<MT> LMIGen(int n){
    int i, j, l, k = n*(n-1)/2+1;
    std::vector<MT> list_Mat;
    MT A;
    list_Mat.push_back(-MT::Identity(n, n));
    for(i = 0; i < n; i++){
        for(j = i+1; j < n; j++){
            A = MT::Zero(n, n);
            A(i,j) = -1;
            A(j,i) = -1;
            list_Mat.push_back(A);
        }
    }
    return list_Mat;
}

template<typename NT, typename Point>
Spectrahedron<Point> prepare_input(int n){
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>     MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 
    int d = n*(n-1)/2;
    Point p(d);
    std::vector<MT> lmi_mat = LMIGen<MT>(n);
    LMI<NT, MT, VT> lmi(lmi_mat);
    Spectrahedron<Point> spectra(lmi);
    spectra.set_interior_point(p);
    return spectra;
}

template <typename NT, typename RNGType, typename Point, typename PointList>
void direct_sampling2(unsigned int n, unsigned int num_points, unsigned int walk_len, PointList &randPoints){

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 

    Spectrahedron<Point> spectra = prepare_input<NT, Point>(n);
    int d = spectra.dimension();
    Point p(d);
    RNGType rng(d);
    std::pair<Point, NT> inner_ball = spectra.ComputeInnerBall();
    NT diameter = 6 * d * inner_ball.second;
    
    for (unsigned int i = 0; i < num_points; ++i){
        p = BilliardWalkSpectra(spectra, p, walk_len, rng, diameter);
        randPoints.push_back(p);
    }
}

template <typename NT, typename WalkType, typename RNGType, typename Point, typename PointList>
void direct_uniform_sampling(int n, int num_points, int walkL, PointList &randPoints, int nburns){

    Spectrahedron<Point> spectra = prepare_input<NT, Point>(n);
    int d = spectra.dimension();
    Point p(d);
    RNGType rng(d);                

    uniform_sampling<WalkType>(randPoints, spectra, rng, walkL, num_points, p, nburns);
}

template <typename NT, typename WalkType, typename RNGType, typename Point, typename PointList>
void direct_gaussian_sampling(int n, int num_points, int walkL, PointList &randPoints, int nburns){

    Spectrahedron<Point> spectra = prepare_input<NT, Point>(n);
    int d = spectra.dimension();
    Point p(d);
    RNGType rng(d);
    NT variance = 1.0;
    NT a = NT(1) / (NT(2) * variance);

    gaussian_sampling<WalkType>(randPoints, spectra, rng, walkL, num_points, a, p, nburns);
}

// void sample_using_gaussian_billiard_walk(HPOLYTOPE const& HP, RNGType& rng, unsigned int walk_len, unsigned int num_points) {
//     std::string filename = "gaussian_billiard_walk.txt";
//     typedef MultivariateGaussianRandomPointGenerator<GaussianAcceleratedWalkType> Generator;

//     std::vector<Point> randPoints;
//     Point q(HP.dimension()); // origin

//     // ----------- Get inscribed ellipsoid --------------------------------
//     typedef Ellipsoid<Point> EllipsoidType;
//     unsigned int max_iter = 150;
//     NT tol = std::pow(10, -6.0), reg = std::pow(10, -4.0);
//     VT x0 = q.getCoefficients();
//     std::pair<std::pair<MT, VT>, bool> inscribed_ellipsoid_res = max_inscribed_ellipsoid<MT>(HP.get_mat(),
//                                                                                              HP.get_vec(),
//                                                                                              x0,
//                                                                                              max_iter,
//                                                                                              tol,
//                                                                                              reg);
//     if (!inscribed_ellipsoid_res.second) // not converged
//         throw std::runtime_error("max_inscribed_ellipsoid not converged");

//     MT A_ell = inscribed_ellipsoid_res.first.first.inverse();
//     EllipsoidType inscribed_ellipsoid(A_ell);
//     // --------------------------------------------------------------------

//     Generator::apply(HP, q, inscribed_ellipsoid, num_points, walk_len,
//                      randPoints, push_back_policy, rng);
//     write_to_file(filename, randPoints);
// }


template <typename NT, typename WalkType, typename RNGType>
void naive_test(int n, int num_points, int walkL, int nreflex){
    typedef Cartesian<NT>                                               Kernel;
    typedef typename Kernel::Point                                      Point;
    typedef std::vector<Point>                                          PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 
    typedef typename WalkType::template Walk<Spectrahedron<Point>, RNGType>   RandomWalk;
    typedef RandomPointGenerator<RandomWalk>                            Generator;

    PointList randPoints;
    
    std::cout << "Direct implementation : " << n << std::endl;
    
    auto start = std::chrono::steady_clock::now();

    direct_uniform_sampling<NT, WalkType, RNGType, Point>(n, num_points, walkL, randPoints, 0);
    
    auto end = std::chrono::steady_clock::now();

    double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << "(ms)" << std::endl;
    
    PointList randPoints2;

    std::cout << "Direct implementation 2 : ";
    
    start = std::chrono::steady_clock::now();
    
    direct_sampling2<NT, RNGType, Point, PointList>(n, num_points, walkL, randPoints2);
    

    end = std::chrono::steady_clock::now();

    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << "(ms)" << std::endl;

    // write_to_file("sampling.txt", randPoints);
    // check_output<VT,MT,PointList>(randPoints, num_points, n);
}
