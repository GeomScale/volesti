#ifndef COMMON_HPP
#define COMMON_HPP
#include "diagnostics/diagnostics.hpp"
#include <unsupported/Eigen/SparseExtra>

inline bool exists_check(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}
/*Problem on the form X=[A|b] bounds=[lb|ub] */
template< typename SpMat, typename VT>
void load_crhmc_problem(SpMat &A, VT &b, VT &lb, VT &ub, int &dimension,
                        std::string problem_name) {
   {
    std::string fileName("./data/");
    fileName.append(problem_name);
    fileName.append(".mm");
    if(!exists_check(fileName)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat X;
    loadMarket(X, fileName);
    int m = X.rows();
    dimension = X.cols() - 1;
    A = X.leftCols(dimension);
    b = VT(X.col(dimension));
  }
  {
    std::string fileName("./data/");
    fileName.append(problem_name);
    fileName.append("_bounds.mm");
    if(!exists_check(fileName)){
      std::cerr<<"Problem does not exist.\n";
      exit(1);}
    SpMat bounds;
    loadMarket(bounds, fileName);
    lb = VT(bounds.col(0));
    ub = VT(bounds.col(1));
  }
}
template <typename NT, typename VT, typename MT>
NT max_interval_psrf(MT &samples) {
  NT max_psrf = NT(0);
  VT intv_psrf = interval_psrf<VT, NT, MT>(samples);
  unsigned int d = intv_psrf.rows();
  for (unsigned int i = 0; i < d; i++) {
    if (intv_psrf(i) > max_psrf)
      max_psrf = intv_psrf(i);
  }
  return max_psrf;
}
template <typename MT, typename VT, typename NT, typename StreamType>
void diagnose(MT &samples, StreamType &stream) {
  unsigned int min_ess = 0;
  print_diagnostics<NT, VT, MT>(samples, min_ess, stream);
  NT max_psrf = max_interval_psrf<NT, VT, MT>(samples);
  stream << "max_psrf: " << max_psrf << std::endl;
  stream << "min ess " << min_ess << std::endl;
}

#endif
