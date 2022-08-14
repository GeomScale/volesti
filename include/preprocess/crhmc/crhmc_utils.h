// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file
#ifndef CRHMC_UTILS_H
#define CRHMC_UTILS_H
#include "Eigen/Eigen"
#include <algorithm>
#include <vector>

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

template <typename SparseMatrixType, typename Type>
void remove_zero_rows(SparseMatrixType &A) {
  std::vector<Eigen::Triplet<Type>> tripletList;
  unsigned Ndata = A.cols();
  unsigned Nbins = A.rows();
  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename SparseMatrixType::InnerIterator it(A, k); it; ++it) {
      tripletList.push_back(
          Eigen::Triplet<Type>(it.row(), it.col(), it.value()));
    }
  }
  // get which rows are empty
  std::vector<bool> has_value(Nbins, false);
  for (auto tr : tripletList)
    has_value[tr.row()] = true;

  if (std::all_of(has_value.begin(), has_value.end(),
                  [](bool v) { return v; })) {
    return;
  }
  // create map from old to new indices
  std::map<unsigned, unsigned> row_map;
  unsigned new_idx = 0;
  for (unsigned old_idx = 0; old_idx < Nbins; old_idx++)
    if (has_value[old_idx])
      row_map[old_idx] = new_idx++;

  // make new triplet list, dropping empty rows
  std::vector<Eigen::Triplet<Type>> newTripletList;
  newTripletList.reserve(Ndata);
  for (auto tr : tripletList)
    newTripletList.push_back(
        Eigen::Triplet<Type>(row_map[tr.row()], tr.col(), tr.value()));

  // form new matrix and return
  SparseMatrixType ret(new_idx, Ndata);
  ret.setFromTriplets(newTripletList.begin(), newTripletList.end());
  A = SparseMatrixType(ret);
}

template <typename SparseMatrixType, typename Type>
void remove_rows(SparseMatrixType &A, std::vector<int> indices) {
  std::vector<Eigen::Triplet<Type>> tripletList;
  unsigned Ndata = A.cols();
  unsigned Nbins = A.rows();
  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename SparseMatrixType::InnerIterator it(A, k); it; ++it) {
      tripletList.push_back(
          Eigen::Triplet<Type>(it.row(), it.col(), it.value()));
    }
  }

  std::vector<bool> notRemoved(Nbins, false);
  for (auto tr : indices)
    notRemoved[tr] = true;

  if (std::all_of(notRemoved.begin(), notRemoved.end(),
                  [](bool v) { return v; })) {
    return;
  }
  // create map from old to new indices
  std::map<unsigned, unsigned> row_map;
  unsigned new_idx = 0;
  for (unsigned old_idx = 0; old_idx < Nbins; old_idx++)
    if (notRemoved[old_idx])
      row_map[old_idx] = new_idx++;

  // make new triplet list, dropping empty rows
  std::vector<Eigen::Triplet<Type>> newTripletList;
  newTripletList.reserve(Ndata);
  for (auto tr : tripletList)
    newTripletList.push_back(
        Eigen::Triplet<Type>(row_map[tr.row()], tr.col(), tr.value()));

  // form new matrix and return
  SparseMatrixType ret(new_idx, Ndata);
  ret.setFromTriplets(newTripletList.begin(), newTripletList.end());
  A = SpMat(ret);
}

#endif
