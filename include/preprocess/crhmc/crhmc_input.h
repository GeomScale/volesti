#ifndef CRHMC_INPUT_H
#define CRHMC_INPUT_H
#include "Eigen/Eigen"
#include "opts.h"

template <typename MatrixType, typename Type> class crhmc_input {
  typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VT;

public:
  MatrixType Aineq;
  VT bineq;
  MatrixType Aeq;
  VT beq;
  opts<Type> options;
  VT lb;
  VT ub;
  crhmc_input(int dimension) {
    Aineq.resize(0, dimension);
    Aeq.resize(0, dimension);
    bineq.resize(0, 1);
    beq.resize(0, 1);
    lb = -VT::Ones(dimension) * std::numeric_limits<Type>::infinity();
    ub = VT::Ones(dimension) * std::numeric_limits<Type>::infinity();
  }
};
#endif
