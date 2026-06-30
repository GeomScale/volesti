// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <list>
#include <stdexcept>

// volesti headers
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "random.hpp"
#include "random/uniform_real_distribution.hpp"
#include "random_walks/random_walks.hpp"
#include "sampling/sampling.hpp"

using namespace Rcpp;

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> SpectrahedronType;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

// [[Rcpp::export]]
SEXP sample_spectrahedra_rcpp(
    SEXP spectrahedra_ptr,
    int n,
    const std::string& walk_type = "RDHR",
    int walk_length = 10,
    int seed = 1,
    Rcpp::Nullable<Rcpp::NumericVector> starting_point = R_NilValue,
    int n_burns = 0) {
  
  try {
    if (n <= 0) {
      throw std::invalid_argument("n must be positive");
    }
    if (walk_length < 1) {
      throw std::invalid_argument("walk_length must be >= 1");
    }
    if (n_burns < 0) {
      throw std::invalid_argument("n_burns cannot be negative");
    }
    
    Rcpp::XPtr<SpectrahedronType> S_ptr(spectrahedra_ptr);
    if (!S_ptr) {
      throw std::runtime_error("Invalid spectrahedron pointer");
    }
    
    SpectrahedronType& S = *S_ptr;
    int dim = S.dimension();
    
    RNGType rng(dim, seed);
    
    Point start_point(dim);
    if (starting_point.isNotNull()) {
      Rcpp::NumericVector sp = Rcpp::as<Rcpp::NumericVector>(starting_point);
      if (sp.size() != dim) {
        throw std::invalid_argument("starting_point dimension mismatch");
      }
      for (int i = 0; i < dim; ++i) {
        start_point.set_coord(i, sp[i]);
      }
    } else {
      for (int i = 0; i < dim; ++i) {
        start_point.set_coord(i, 0.0);
      }
    }
    
    std::list<Point> samples;
    
    if (walk_type == "RDHR") {
      typedef RandomDirichletHitAndRunWalk WalkType;
      uniform_sampling<WalkType>(
        samples, S, rng, walk_length, n, start_point, n_burns
      );
      
    } else if (walk_type == "CDHR") {
      typedef CoordinateDirectedHitAndRunWalk WalkType;
      uniform_sampling<WalkType>(
        samples, S, rng, walk_length, n, start_point, n_burns
      );
      
    } else if (walk_type == "BILLIARD") {
      typedef AcceleratedBilliardWalk WalkType;
      uniform_sampling<WalkType>(
        samples, S, rng, walk_length, n, start_point, n_burns
      );
      
    } else {
      throw std::invalid_argument("Unknown walk type: " + walk_type);
    }
    
    int n_samples = samples.size();
    if (n_samples == 0) {
      throw std::runtime_error("No samples were generated");
    }
    
    MT sample_matrix(dim, n_samples);
    int col = 0;
    for (const auto& p : samples) {
      sample_matrix.col(col) = p.getCoefficients();
      col++;
    }
    
    return Rcpp::wrap(sample_matrix);
    
  } catch(const std::exception& e) {
    std::string error_msg = "C++ error: ";
    error_msg += e.what();
    stop(error_msg);
    return R_NilValue;
    
  } catch(...) {
    stop("Unknown C++ exception");
    return R_NilValue;
  }
}
