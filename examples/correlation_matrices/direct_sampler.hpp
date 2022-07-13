// Direct implementation: call to sampling algorithms for spectrahedra in volesti
#include "sampling/sphere.hpp"
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

template <typename NT, typename Point, typename RNGType>
Point BilliardWalkSpectra(Spectrahedron<Point> &P, Point& q, unsigned int const& walk_length, unsigned int nreflex, RNGType &rng, NT const _Len){
    unsigned int k = P.dimension();
    NT L, tau;
    Point p = q, v;
    int flag = -1;
    for (unsigned int j=0; j<walk_length; ++j){
        L = rng.sample_urdist() * _Len;
        v = GetDirection<Point>::apply(k, rng);
        Point p0 = p;
        int it = 0;
        while (it < nreflex)
        {
            tau = P.positiveLinearIntersection(p.getCoefficients(), v.getCoefficients());
            if (L <= tau) {
                p += (L * v);
                break;
            }
            tau = 0.995 * tau; // 0.995: to approximate boundary points?
            p += tau * v; // A point (almost) on the boundary
            L -= tau;
            P.compute_reflection(v, p, flag); // reflection(P, p, v, pbpair.second);
            it++;
        }
        if (it == nreflex){
            p = p0;
        }
    }
    return p;
}

template <typename NT, typename RNGType, typename Point, typename PointList>
void direct_sampling2(unsigned int n, unsigned int num_points, unsigned int walk_len, unsigned int nreflex, PointList &randPoints){

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 

    int d = n*(n-1)/2; // Dimension: k = n(n-1)/2
    RNGType rng(d);
    Point p(d); // Initial interior point (origin: identity matrix)

    // Create the spectrahedron
    std::vector<MT> lmi_mat = LMIGen<MT>(n);
    LMI<NT, MT, VT> lmi(lmi_mat);
    Spectrahedron<Point> spectra(lmi);
    spectra.set_interior_point(p);
    std::pair<Point, NT> inner_ball = spectra.ComputeInnerBall();
    NT diameter = 6 * d * inner_ball.second;
    
    for (unsigned int i = 0; i < num_points; ++i){
        p = BilliardWalkSpectra(spectra, p, walk_len, nreflex, rng, diameter);
        randPoints.push_back(p);
    }
}

template <typename NT, typename WalkType, typename RNGType, typename Point, typename PointList>
void direct_uniform_sampling(int n, int num_points, int walkL, int nreflex, PointList &randPoints){

    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 
    typedef typename WalkType::template Walk<Spectrahedron<Point>, RNGType>   RandomWalk;
    typedef RandomPointGenerator<RandomWalk>                            Generator;

    int d = n*(n-1)/2;
    RNGType rng(d);
    Point p(d); // Initial interior point (origin: identity matrix)

    // Create the spectrahedron
    std::vector<MT> lmi_mat = LMIGen<MT>(n);
    LMI<NT, MT, VT> lmi(lmi_mat);
    Spectrahedron<Point> spectra(lmi);
    spectra.set_interior_point(p);
    
    PushBackWalkPolicy push_back_policy;
    
    Generator::apply(spectra, p, num_points, walkL,
                     randPoints, push_back_policy, rng);
}

template <typename NT, typename WalkType, typename RNGType>
void naive_test(int n, int num_points, int walkL, int nreflex){
    typedef Cartesian<NT>                                               Kernel;
    typedef typename Kernel::Point                                      Point;
    typedef std::vector<Point>                                          PointList;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>             MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                          VT; 

    PointList randPoints;
    
    std::cout << "Direct implementation : ";
    
    auto start = std::chrono::steady_clock::now();

    direct_sampling2<NT, RNGType, Point, PointList>(n, num_points, walkL, nreflex, randPoints);
    // direct_sampling2<NT, RNGType, Point, PointList>(n, num_points, walkL, nreflex, randPoints);

    auto end = std::chrono::steady_clock::now();

    double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << "(ms)" << std::endl;

    PointList randPoints2;

    std::cout << "Direct implementation 2 : ";
    
    start = std::chrono::steady_clock::now();

    direct_uniform_sampling<NT, WalkType, RNGType, Point, PointList>(n, num_points, walkL, nreflex, randPoints2);

    end = std::chrono::steady_clock::now();

    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << time << "(ms)" << std::endl;
    write_to_file("uniform_sampling.txt", randPoints);
    write_to_file("uniform_sampling2.txt", randPoints2);
    // for(int i = 0; i < num_points ; ++i){
    //     if(!membership<VT,MT>(randPoints2[i].getCoefficients(), n)){
    //         std::cout << "ALERT\n";
    //     }else{ std::cout << "OK\n";}
            

    // }
}
