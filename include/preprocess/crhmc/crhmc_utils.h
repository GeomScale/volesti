#ifndef CRHMC_UTILS_H
#define CRHMC_UTILS_H
#include "Eigen/Eigen"

template <typename Func> struct lambda_as_visitor_wrapper : Func {
  lambda_as_visitor_wrapper(const Func &f) : Func(f) {}
  template <typename S, typename I> void init(const S &v, I i, I j) {
    return Func::operator()(v, i, j);
  }
};

template <typename Mat, typename Func>
void visit_lambda(const Mat &m, const Func &f) {
  lambda_as_visitor_wrapper<Func> visitor(f);
  m.visit(visitor);
}

template <typename SparseMatrixType>
void sparse_stack_v(const SparseMatrixType &top, const SparseMatrixType &bottom,
                    SparseMatrixType &stacked) {
  assert(top.cols() == bottom.cols());
  stacked.resize(top.rows() + bottom.rows(), top.cols());
  stacked.resizeNonZeros(top.nonZeros() + bottom.nonZeros());

  int i = 0;

  for (int col = 0; col < top.cols(); col++) {
    stacked.outerIndexPtr()[col] = i;

    for (int j = top.outerIndexPtr()[col]; j < top.outerIndexPtr()[col + 1];
         j++, i++) {
      stacked.innerIndexPtr()[i] = top.innerIndexPtr()[j];
      stacked.valuePtr()[i] = top.valuePtr()[j];
    }

    for (int j = bottom.outerIndexPtr()[col];
         j < bottom.outerIndexPtr()[col + 1]; j++, i++) {
      stacked.innerIndexPtr()[i] = (int)top.rows() + bottom.innerIndexPtr()[j];
      stacked.valuePtr()[i] = bottom.valuePtr()[j];
    }
  }
  stacked.outerIndexPtr()[top.cols()] = i;
}

template <typename SparseMatrixType>
void sparse_stack_h(const SparseMatrixType &left, const SparseMatrixType &right,
                    SparseMatrixType &stacked) {
  assert(left.rows() == right.rows());

  stacked.resize(left.rows(), left.cols() + right.cols());
  stacked.resizeNonZeros(left.nonZeros() + right.nonZeros());

  std::copy(left.innerIndexPtr(), left.innerIndexPtr() + left.nonZeros(),
            stacked.innerIndexPtr());
  std::copy(right.innerIndexPtr(), right.innerIndexPtr() + right.nonZeros(),
            stacked.innerIndexPtr() + left.nonZeros());

  std::copy(left.valuePtr(), left.valuePtr() + left.nonZeros(),
            stacked.valuePtr());
  std::copy(right.valuePtr(), right.valuePtr() + right.nonZeros(),
            stacked.valuePtr() + left.nonZeros());

  std::copy(left.outerIndexPtr(), left.outerIndexPtr() + left.cols(),
            stacked.outerIndexPtr()); // dont need the last entry of
                                      // A.outerIndexPtr() -- total length is
                                      // AB.cols() + 1 = A.cols() + B.cols() + 1
  std::transform(right.outerIndexPtr(),
                 right.outerIndexPtr() + right.cols() + 1,
                 stacked.outerIndexPtr() + left.cols(),
                 [&](int i) { return i + left.nonZeros(); });
}

template <typename SparseMatrixType>
void sparse_stack_h_inplace(SparseMatrixType &left,
                            const SparseMatrixType &right) {
  assert(left.rows() == right.rows());

  const int leftcol = (int)left.cols();
  const int leftnz = (int)left.nonZeros();

  left.conservativeResize(left.rows(), left.cols() + right.cols());
  left.resizeNonZeros(left.nonZeros() + right.nonZeros());

  std::copy(right.innerIndexPtr(), right.innerIndexPtr() + right.nonZeros(),
            left.innerIndexPtr() + leftnz);
  std::copy(right.valuePtr(), right.valuePtr() + right.nonZeros(),
            left.valuePtr() + leftnz);
  std::transform(
      right.outerIndexPtr(), right.outerIndexPtr() + right.cols() + 1,
      left.outerIndexPtr() + leftcol, [&](int i) { return i + leftnz; });
}

#endif
