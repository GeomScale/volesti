// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <stdexcept>

// volesti headers
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "random.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_sequence_of_balls.hpp"

using namespace Rcpp;

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> SpectrahedronType;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

// [[Rcpp::export]]
double volume_spectrahedra_rcpp(
    SEXP spectrahedra_ptr,
    int walk_length = 10,
    double tolerance = 0.1,
    int seed = 1,
    const std::string& volume_method = "CB") {
  
  try {
    if (walk_length < 1) {
      throw std::invalid_argument("walk_length must be >= 1");
    }
    if (tolerance <= 0) {
      throw std::invalid_argument("tolerance must be positive");
    }
    
    Rcpp::XPtr<SpectrahedronType> S_ptr(spectrahedra_ptr);
    
    if (!S_ptr) {
      throw std::runtime_error("Invalid spectrahedron pointer");
    }
    
    SpectrahedronType& S = *S_ptr;
    int dim = S.dimension();
    
    RNGType rng(dim, seed);
    
    Point interior_point(dim);
    for (int i = 0; i < dim; ++i) {
      interior_point.set_coord(i, 0.0);
    }
    
    NT volume = 0.0;
    
    if (volume_method == "CB") {
      typedef AcceleratedBilliardWalk WalkType;
      
      auto result = volume_cooling_balls<WalkType, RNGType>(
        S, interior_point, walk_length, tolerance
      );
      
      volume = result.second;
      
    } else if (volume_method == "SOB") {
      typedef AcceleratedBilliardWalk WalkType;
      
      auto result = volume_sequence_of_balls<WalkType, RNGType>(
        S, interior_point, walk_length, tolerance
      );
      
      volume = result.second;
      
    } else {
      throw std::invalid_argument("Unknown volume method: " + volume_method);
    }
    
    return volume;
    
  } catch(const std::exception& e) {
    std::string error_msg = "C++ error: ";
    error_msg += e.what();
    stop(error_msg);
    return NA_REAL;
    
  } catch(...) {
    stop("Unknown C++ exception");
    return NA_REAL;
  }
}
