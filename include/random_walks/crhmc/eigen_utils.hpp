#ifndef EIGEN_UTILS_HPP
#define EIGEN_UTILS_HPP
#include "Eigen/Eigen"

using VT =Eigen::Array<double,Eigen::Dynamic,1>;
using MT =Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic>;

inline MT operator+(const MT &M, const VT v) {
  return M.colwise()+v;
}
inline MT operator-(const MT &M, const VT v) {
  return M.colwise()-v;
}
inline MT operator-( const VT v,const MT &M) {
  return (-M).colwise()+v;
}

#endif
