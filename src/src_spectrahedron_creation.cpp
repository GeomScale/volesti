// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <fstream>
#include <stdexcept>
#include <vector>

// volesti headers
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "convex_bodies/spectrahedra/LMI.h"
#include "SDPAFormatManager.h"

using namespace Rcpp;

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron<Point> SpectrahedronType;
typedef LMI<NT, MT, VT> LMIType;

// [[Rcpp::export]]
SEXP load_spectrahedron_from_sdpa_cpp(const std::string& filename) {
  try {
    std::ifstream test_file(filename);
    if (!test_file.is_open()) {
      throw std::runtime_error("Cannot open file: " + filename);
    }
    test_file.close();
    
    auto* S = new SpectrahedronType();
    Point objFunction;
    
    SdpaFormatManager<NT> manager;
    std::ifstream infile(filename);
    
    manager.loadSDPAFormatFile(infile, *S, objFunction);
    
    infile.close();
    
    Rcpp::XPtr<SpectrahedronType> ptr(S, true);
    
    return ptr;
    
  } catch(const std::exception& e) {
    std::string error_msg = "Error loading SDPA file: ";
    error_msg += e.what();
    stop(error_msg);
    return R_NilValue;
  } catch(...) {
    stop("Unknown C++ exception");
    return R_NilValue;
  }
}

// [[Rcpp::export]]
int get_spectrahedron_dimension_cpp(SEXP ptr_sexp) {
  try {
    Rcpp::XPtr<SpectrahedronType> ptr(ptr_sexp);
    if (!ptr) {
      throw std::runtime_error("Invalid pointer");
    }
    return ptr->dimension();
    
  } catch(const std::exception& e) {
    std::string error_msg = "Error getting dimension: ";
    error_msg += e.what();
    stop(error_msg);
    return -1;
  }
}

// [[Rcpp::export]]
SEXP create_spectrahedron_from_matrices_cpp(
    Rcpp::NumericMatrix L0_r,
    Rcpp::List L_list_r) {
  
  try {
    MT L0 = Rcpp::as<MT>(L0_r);
    
    if (L0.rows() != L0.cols()) {
      throw std::invalid_argument("L0 must be square");
    }
    
    int m = L0.rows();
    int n = L_list_r.size();
    
    std::vector<MT> L_matrices;
    for (int i = 0; i < n; ++i) {
      MT L_i = Rcpp::as<MT>(L_list_r[i]);
      if (L_i.rows() != m || L_i.cols() != m) {
        throw std::invalid_argument("All matrices must have same dimension");
      }
      L_matrices.push_back(L_i);
    }
    
    auto* S = new SpectrahedronType();
    
    Rcpp::XPtr<SpectrahedronType> ptr(S, true);
    return ptr;
    
  } catch(const std::exception& e) {
    std::string error_msg = "Error creating spectrahedron: ";
    error_msg += e.what();
    stop(error_msg);
    return R_NilValue;
  }
}

// [[Rcpp::export]]
bool is_point_in_spectrahedron_cpp(SEXP ptr_sexp, Rcpp::NumericVector point_r) {
  try {
    Rcpp::XPtr<SpectrahedronType> ptr(ptr_sexp);
    if (!ptr) {
      throw std::runtime_error("Invalid pointer");
    }
    
    SpectrahedronType& S = *ptr;
    
    Point p(S.dimension());
    if (point_r.size() != S.dimension()) {
      throw std::invalid_argument("Point dimension mismatch");
    }
    
    for (int i = 0; i < S.dimension(); ++i) {
      p.set_coord(i, point_r[i]);
    }
    
    return (S.is_in(p) == 0);
    
  } catch(const std::exception& e) {
    std::string error_msg = "Error checking point: ";
    error_msg += e.what();
    stop(error_msg);
    return false;
  }
}
