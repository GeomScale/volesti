// A lot of code here is not my original idea. I borrow the structure and 
// a lot of ideas from the source code of volesti

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <boost/random.hpp>
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/spectrahedra/spectrahedron.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VT;

typedef Cartesian<double> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> spectrahedron;
typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;

// Assume that A is already symmetric
bool isPosSemidefinite(MT A){
    // Eigen::LLT<MT> A_llt(A);
    // if(A_llt.info() != Eigen::NumericalIssue) return true;
    // return false;
    Eigen::LDLT<MT> A_ldlt(A);
    if (A_ldlt.info() != Eigen::NumericalIssue && A_ldlt.isPositive())
        return true;
    return false;
}

// A small function to build a correlation matrix from a vector of entries
MT rebuildMatrix(const Point &p, const unsigned int n){
    VT xvector = p.getCoefficients();
    MT A = MT::Identity(n,n);
    double coeff;
    for(int i = 0; i < n ; ++i){
        for(int j = i+1; j < n; ++j){
            int ind = ((((n<<1)-i-2)*(i+1)) >> 1)  + j - n;
            coeff = xvector[ind];
            A(i,j) = coeff;
            A(j,i) = coeff;
        }
    }
    return A;
}

std::vector<MT> myLMIGenerator(int n){
    int i, j, l, k = n*(n-1)/2+1;
    std::vector<MT> list_Mat;
    MT A;
    list_Mat.push_back(MT::Identity(n, n));
    for(i = 0; i < n; i++){
        for(j = i+1; j < n; j++){
            A = MT::Zero(n, n);
            A(i,j) = 1;
            A(j,i) = 1;
            list_Mat.push_back(A);
        }
    }
    return list_Mat;
}

// Auxiliary geometric functions
Point getDirection(unsigned int const& dim, RNGType &rng, bool normalize=true){
    double normal = 0.;
    Point p(dim);
    double* data = p.pointerToData();

    for (unsigned int i=0; i<dim; ++i){
        *data = rng.sample_ndist();
        normal += *data * *data;
        data++;
    }

    normal = 1./std::sqrt(normal);
    if (normalize) p *= normal;
    return p;
}

bool membership(spectrahedron &spectra, const VT &xvector, const unsigned int n, const unsigned int k){
    // std::vector<double>::iterator it = xvector.begin();
    for(int i = 0; i < k; ++i)
        if((xvector(i) > 1) || (xvector(i) < -1)) return false;
    MT A = rebuildMatrix(xvector, n);
    if(isPosSemidefinite(A)) return true;
    return false;
}

std::pair<double, int> intersection(spectrahedron &P, const Point &x, const Point &v, const unsigned int k){
    double tau, tmp;
    int j = 0;
    if(v[0] > 0){
        tau = (1-x[0])/v[0];   
    }else{
        tau = -(1 + x[0])/v[0];
    }
    for(int i = 1; i < k; ++i){
        if(v[i] > 0){
            tmp = (1 - x[i])/v[i];
        }else{
            tmp = -(1 + x[i])/v[i];
        }
        if(tau > tmp){
            tau = tmp;
            j = i;
        }
    }
    tmp = P.positiveLinearIntersection(x.getCoefficients(), v.getCoefficients());
    if(tau > tmp){
         tau = tmp;
        j = -1;
    }
    std::pair<double, int> res(tau,j);
    return res;
}

void reflection(spectrahedron P, Point &p, Point &v, const int flag){
    if(flag != -1){
        v.set_coord(flag, - v.getCoefficients()(flag));
        return;
    }
    P.compute_reflection(v, p, flag);
}

// This function is taken and simplified from uniform_billiard_walk.hpp

Point BilliardWalkSpectra(spectrahedron &P, Point& q, unsigned int const& walk_length, unsigned int nreflex, RNGType &rng, double const _Len){
    unsigned int k = P.dimension();
    double L, tau;
    Point p = q, v;
    std::pair<double, int> pbpair;
    for (unsigned int j=0; j<walk_length; ++j){
        L = rng.sample_urdist() * _Len;
        v = getDirection(k, rng);
        Point p0 = p;
        int it = 0;
        while (it < nreflex)
        {
            pbpair = intersection(P, p, v, k);
            tau = pbpair.first;
            if (L <= tau) {
                p += (L * v);
                break;
            }
            tau = 0.995 * tau; // 0.995: to approximate boundary points?
            p += tau * v; // A point (almost) on the boundary
            L -= tau;
            reflection(P, p, v, pbpair.second);
            it++;
        }
        if (it == nreflex){
            p = p0;
        }
    }
    return p;
}

std::vector<Point> uniform_correl_matrix_sampling(unsigned int n, unsigned int num_points, unsigned int walk_len=10, unsigned int nreflex = 10){
    int k = n*(n-1)/2; // Dimension: k = n(n-1)/2
    RNGType rng(k);
    Point p(k); // Initial interior point (origin: identity matrix)
    std::vector<Point> randPoints;

    // Create the spectrahedron
    std::vector<MT> lmi_mat = myLMIGenerator(n);
    LMI<double, MT, VT> lmi(lmi_mat);
    spectrahedron spectra(lmi);
    spectra.set_interior_point(p);
    std::pair<Point, double> inner_ball = spectra.ComputeInnerBall();
    double diameter = 6 * k * inner_ball.second;
    
    for (unsigned int i = 0; i < num_points; ++i){
        p = BilliardWalkSpectra(spectra, p, walk_len, nreflex, rng, diameter);
        randPoints.push_back(p);
    }
    return randPoints;
}

// Test 3:

void Test(unsigned int num_points=1000, unsigned int walk_len=10, unsigned int nreflex = 10){
    int n = 8;
    double time;
    while(n < 9){

        auto start = std::chrono::steady_clock::now();

        uniform_correl_matrix_sampling(n, num_points, walk_len, nreflex);

        auto end = std::chrono::steady_clock::now();
        time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << n << " : " << time << "(ms)" << std::endl;
        ++n;
    }
}

int main(){
    srand((unsigned) time(NULL));
    int n = 8;
    Test();
    return 0;
}