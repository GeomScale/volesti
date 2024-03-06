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
#include <unsupported/Eigen/SparseExtra>
#include "PackedCSparse/SparseMatrix.h"
#include <algorithm>
#include <vector>

template <typename Func>
struct lambda_as_visitor_wrapper : Func
{
  lambda_as_visitor_wrapper(const Func &f) : Func(f) {}
  template <typename S, typename I>
  void init(const S &v, I i, I j)
  {
    return Func::operator()(v, i, j);
  }
};

template <typename Mat, typename Func>
void visit_lambda(const Mat &m, const Func &f)
{
  lambda_as_visitor_wrapper<Func> visitor(f);
  m.visit(visitor);
}

template <typename SparseMatrixType>
void sparse_stack_v(const SparseMatrixType &top, const SparseMatrixType &bottom,
                    SparseMatrixType &stacked)
{
  assert(top.cols() == bottom.cols());
  stacked.resize(top.rows() + bottom.rows(), top.cols());
  stacked.resizeNonZeros(top.nonZeros() + bottom.nonZeros());

  int i = 0;

  for (int col = 0; col < top.cols(); col++)
  {
    stacked.outerIndexPtr()[col] = i;

    for (int j = top.outerIndexPtr()[col]; j < top.outerIndexPtr()[col + 1];
         j++, i++)
    {
      stacked.innerIndexPtr()[i] = top.innerIndexPtr()[j];
      stacked.valuePtr()[i] = top.valuePtr()[j];
    }

    for (int j = bottom.outerIndexPtr()[col];
         j < bottom.outerIndexPtr()[col + 1]; j++, i++)
    {
      stacked.innerIndexPtr()[i] = (int)top.rows() + bottom.innerIndexPtr()[j];
      stacked.valuePtr()[i] = bottom.valuePtr()[j];
    }
  }
  stacked.outerIndexPtr()[top.cols()] = i;
}

template <typename SparseMatrixType>
void sparse_stack_h(const SparseMatrixType &left, const SparseMatrixType &right,
                    SparseMatrixType &stacked)
{
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
                 [&](int i)
                 { return i + left.nonZeros(); });
}
#include <unsupported/Eigen/SparseExtra>

template <typename SparseMatrixType>
void sparse_stack_h_inplace(SparseMatrixType &left,
                            const SparseMatrixType &right)
{
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
      left.outerIndexPtr() + leftcol, [&](int i)
      { return i + leftnz; });
}

template <typename SparseMatrixType, typename VectorType, typename Type>
void remove_zero_rows(SparseMatrixType &A, VectorType &b) {
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
    if (has_value[old_idx]) {
      row_map[old_idx] = new_idx;
      b(new_idx) = b(old_idx);
      new_idx++;
    }
  // make new triplet list, dropping empty rows
  std::vector<Eigen::Triplet<Type>> newTripletList;
  newTripletList.reserve(Ndata);
  for (auto tr : tripletList)
    newTripletList.push_back(
        Eigen::Triplet<Type>(row_map[tr.row()], tr.col(), tr.value()));

  // form new matrix and return
  A.resize(new_idx, Ndata);
  A.setFromTriplets(newTripletList.begin(), newTripletList.end());
  b.conservativeResize(new_idx);
}

template <typename SparseMatrixType, typename Type>
void remove_rows(SparseMatrixType &A, std::vector<bool> &notRemoved) {
  unsigned Ndata = A.cols();
  unsigned Nbins = A.rows();
  // create map from old to new indices
  std::map<unsigned, unsigned> row_map;
  unsigned new_idx = 0;
  for (unsigned old_idx = 0; old_idx < Nbins; old_idx++)
    if (notRemoved[old_idx])
      row_map[old_idx] = new_idx++;

  std::vector<Eigen::Triplet<Type>> tripletList;
  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename SparseMatrixType::InnerIterator it(A, k); it; ++it) {
      if (notRemoved[it.row()]) {
        tripletList.push_back(
            Eigen::Triplet<Type>(row_map[it.row()], it.col(), it.value()));
      }
    }
  }
  // form new matrix and return
  A.resize(new_idx, Ndata);
  A.setFromTriplets(tripletList.begin(), tripletList.end());
}

template <typename SparseMatrixType, typename VectorType, typename Type>
std::pair<VectorType, VectorType> colwiseMinMax(SparseMatrixType const &A)
{
  int n = A.cols();
  VectorType cmax(n);
  VectorType cmin(n);
  for (int k = 0; k < A.outerSize(); ++k)
  {
    Type minv = +std::numeric_limits<Type>::max();
    Type maxv = std::numeric_limits<Type>::lowest();
    for (typename SparseMatrixType::InnerIterator it(A, k); it; ++it)
    {
      minv = std::min(minv, it.value());
      maxv = std::max(maxv, it.value());
    }
    cmin(k) = minv;
    cmax(k) = maxv;
  }
  return std::make_pair(cmin, cmax);
}
template <typename VectorType>
void nextpow2(VectorType &a)
{
  a = (a.array() == 0).select(1, a);
  a = (((a.array().log()) / std::log(2)).ceil()).matrix();
  a = pow(2, a.array()).matrix();
}
template <typename SparseMatrixType, typename VectorType, typename Type>
std::pair<VectorType, VectorType> gmscale(SparseMatrixType &Asp,
                                          const Type scltol)
{
  using Diagonal_MT = Eigen::DiagonalMatrix<Type, Eigen::Dynamic>;
  int m = Asp.rows();
  int n = Asp.cols();
  SparseMatrixType A = Asp.cwiseAbs();
  A.makeCompressed();
  int maxpass = 10;
  Type aratio = 1e+50;
  Type sratio;
  Type damp = 1e-4;
  Type small = 1e-8;
  VectorType rscale = VectorType ::Ones(m, 1);
  VectorType cscale = VectorType ::Ones(n, 1);
  VectorType cmax;
  VectorType cmin;
  VectorType rmax;
  VectorType rmin;
  VectorType eps = VectorType ::Ones(n, 1) * 1e-12;
  SparseMatrixType SA;
  for (int npass = 0; npass < maxpass; npass++)
  {

    rscale = (rscale.array() == 0).select(1, rscale);
    Diagonal_MT Rinv = (rscale.cwiseInverse()).asDiagonal();
    SA = Rinv * A;
    std::tie(cmin, cmax) =
        colwiseMinMax<SparseMatrixType, VectorType, Type>(SA);

    // cmin = (cmin + eps).cwiseInverse();
    sratio = (cmax.cwiseQuotient(cmin)).maxCoeff();

    if (npass > 0)
    {
      cscale = ((cmin.cwiseMax(damp * cmax)).cwiseProduct(cmax)).cwiseSqrt();
    }

    if (npass >= 2 && sratio >= aratio * scltol)
    {
      break;
    }
    aratio = sratio;
    nextpow2(cscale);
    Diagonal_MT Cinv = (cscale.cwiseInverse()).asDiagonal();
    SA = A * Cinv;
    std::tie(rmin, rmax) =
        colwiseMinMax<SparseMatrixType, VectorType, Type>(SA.transpose());
    // rmin = (rmin + eps).cwiseInverse();
    rscale = ((rmin.cwiseMax(damp * rmax)).cwiseProduct(rmax)).cwiseSqrt();
    nextpow2(rscale);
  }
  rscale = (rscale.array() == 0).select(1, rscale);
  Diagonal_MT Rinv = (rscale.cwiseInverse()).asDiagonal();
  SA = Rinv * A;
  std::tie(std::ignore, cscale) =
      colwiseMinMax<SparseMatrixType, VectorType, Type>(SA);
  nextpow2(cscale);
  return std::make_pair(cscale, rscale);
}
template <typename Type>
int doubleVectorEqualityComparison(
    const Type a, const Type b,
    const Type tol = std::numeric_limits<Type>::epsilon())
{
  return (abs(a - b) < tol * (1 + abs(a) + abs(b)));
}

template <typename SparseMatrixType>
std::pair<std::vector<int>, std::vector<int>>
nnzPerColumn(SparseMatrixType const &A, const int threashold)
{
  int n = A.cols();
  std::vector<int> colCounts(n);
  std::vector<int> badCols;
  for (int k = 0; k < A.outerSize(); ++k)
  {
    int nnz = 0;
    for (typename SparseMatrixType::InnerIterator it(A, k); it; ++it)
    {
      if (it.value() != 0)
      {
        nnz++;
      }
    }
    colCounts[k] = nnz;
    if (nnz > threashold)
    {
      badCols.push_back(k);
    }
  }
  return std::make_pair(colCounts, badCols);
}
using PM = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;
template <typename SparseMatrixType>
PM permuteMatAMD(SparseMatrixType const &A)
{
  Eigen::AMDOrdering<int> ordering;
  PM perm;
  ordering(A, perm);
  return perm;
}
template <typename SparseMatrixType>
PM postOrderPerm(SparseMatrixType const &A)
{
  using IndexVector = Eigen::Matrix<int, Eigen::Dynamic, 1>;
  int n = A.rows();
  IndexVector m_etree;
  IndexVector firstRowElt;
  Eigen::internal::coletree(A, m_etree, firstRowElt);
  IndexVector post;
  Eigen::internal::treePostorder(int(A.cols()), m_etree, post);
  PM post_perm(n);
  for (int i = 0; i < n; i++)
    post_perm.indices()(i) = post(i);
  return post_perm;
}

template <typename SparseMatrixType,typename VectorType>
void fillin_reduce(SparseMatrixType &X,VectorType& b){
  SparseMatrixType I = SparseMatrixType(Eigen::VectorXd::Ones(X.rows()).asDiagonal());
  SparseMatrixType XX = X * X.transpose() + I;
  XX.makeCompressed();
  Eigen::SimplicialLDLT<SparseMatrixType, Eigen::Lower,
                        Eigen::AMDOrdering<int>> cholesky;
  cholesky.analyzePattern(XX);
  X = cholesky.permutationP() * X;
  b = cholesky.permutationP() *b;
}
template<typename SparseMatrixType,typename Type,typename IndexType>
PackedCSparse::SparseMatrix<Type,IndexType> transform_format(SparseMatrixType const &mat) {
  PackedCSparse::SparseMatrix<Type, IndexType> A = PackedCSparse::SparseMatrix<Type, IndexType>(mat.rows(), mat.cols(), mat.nonZeros());
  IndexType nnz = 0;
  for (IndexType outeindex = 0; outeindex < mat.outerSize(); ++outeindex) {
    A.p[outeindex] = nnz;
    for (typename SparseMatrixType::InnerIterator it(mat, outeindex); it; ++it) {
      A.i[nnz] = it.row();
      A.x[nnz] = it.value();
      nnz++;
    }
  }
  A.p[A.n] = nnz;
  return A;
}
template<typename MatrixType, typename IndexType>
void copy_indicies(MatrixType& a, MatrixType& b, std::vector<IndexType>const & a_idx, std::vector<IndexType>const & b_idx){
for(int i=0;i<b_idx.size();i++){
  a(a_idx[i])=b(b_idx[i]);
}
}
template<typename MatrixType, typename IndexType>
void copy_indicies(MatrixType& a, MatrixType b, std::vector<IndexType>const & b_idx){
for(int i=0;i<b_idx.size();i++){
  a(i)=b(b_idx[i]);
}
}
template<typename MatrixType, typename IndexType, typename Type>
void set(MatrixType &a, std::vector<IndexType>const & idx, const Type c){
  for(int i=0;i<idx.size();i++){
    a(idx[i])=c;
  }
}
template<typename MatrixType, typename IndexType>
void volesti_saxpy(MatrixType &a,MatrixType const &b,MatrixType const& c, std::vector<IndexType>const & a_idx, std::vector<IndexType>const & b_idx){
for(int i=0;i<b_idx.size();i++){
  a(a_idx[i])=b(b_idx[i])+c(i);
}
}
/*Problem on the form X=[A|b] bounds=[lb|ub] */
template< typename SpMat, typename VT>
void load_problem(SpMat &A, VT &b, VT &lb, VT &ub, int &dimension,
                        std::string problem_name) {
   {
    std::string fileName(problem_name);
    fileName.append(".mm");
    SpMat X;
    loadMarket(X, fileName);
    int m = X.rows();
    dimension = X.cols() - 1;
    A = X.leftCols(dimension);
    b = VT(X.col(dimension));
  }
  {
    std::string fileName(problem_name);
    fileName.append("_bounds.mm");
    SpMat bounds;
    loadMarket(bounds, fileName);
    lb = VT(bounds.col(0));
    ub = VT(bounds.col(1));
  }
}
#endif
