/// Functions to sample correlation matrices w.r.t. a truncated density

#ifndef SAMPLE_CORRELATION_MATRICES_HPP
#define SAMPLE_CORRELATION_MATRICES_HPP

#include <sampling/sampling.hpp>

/// Direct implementation: call to sampling algorithms for spectrahedra in volesti
/// Mainly serves for comparison

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
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1>                  VT;

    int d = n*(n-1)/2;
    Point p(d);
    std::vector<MT> lmi_mat = LMIGen<MT>(n);
    LMI<NT, MT, VT> lmi(lmi_mat);
    Spectrahedron<Point> spectra(lmi);
    spectra.set_interior_point(p);
    return spectra;
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

// New implementations for sampling correlation matrices

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename PointList
>
void uniform_correlation_sampling(const unsigned int &n,
                                    PointList &randPoints,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    unsigned int const& nburns){
    CorreSpectra<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(d);
    RNGType rng(d);
    
    uniform_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, startingPoint,0);
}

template
<
    typename WalkTypePolicy,
    typename PointType,
    typename RNGType,
    typename PointList,
    typename NT
>
void gaussian_correlation_sampling(const unsigned int &n,
                                    PointList &randPoints,
                                    const unsigned int &walkL,
                                    const unsigned int &num_points,
                                    const NT &a,
                                    unsigned int const& nburns = 0){
    CorreSpectra<PointType> P(n);
    const unsigned int d = P.dimension();
    PointType startingPoint(d);
    RNGType rng(d);
    
    gaussian_sampling<WalkTypePolicy>(randPoints, P, rng, walkL, num_points, a, startingPoint, 0);
}

#endif
