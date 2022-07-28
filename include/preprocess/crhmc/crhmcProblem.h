#ifndef CRHMCPROBLEM_H
#define CRHMCPROBLEM_H

#include "Eigen/Eigen"
#include "PackedChol.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "hpolytope.h"
#include "random_walks/random_walks.hpp"
#include "sos/barriers/TwoSidedBarrier.h"
#include <vector>
#include "input_structure.h"

#include "analytic_center.h"
#include "lewis_center.h"
#include "misc.h"
#include "opts.h"
#include "volume_cooling_balls.hpp"
#include "volume_cooling_gaussians.hpp"
#include "volume_sequence_of_balls.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include "crhmc_utils.h"

#ifndef SIMD_LEN
#define SIMD_LEN 0
#endif
const size_t chol_k = (SIMD_LEN == 0) ? 1 : SIMD_LEN;

template <typename Point> class crhmcProblem {
public:
  typedef double NT;
  typedef Cartesian<NT> Kernel;
  typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
  typedef HPolytope<Point> HPOLYTOPE;
  typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
  typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
  typedef Eigen::SparseMatrix<NT> SpMat;
  typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> PM;
  typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IndexVector;
  typedef PackedChol<chol_k, int> CholObj;
  typedef Eigen::Triplet<double> Triple;
  typedef TwoSidedBarrier<NT> SimpleBarrier;
  using Tx2 = FloatArray<double, chol_k>;
  typedef opts<NT> Opts;
  typedef Eigen::DiagonalMatrix<NT, Eigen::Dynamic> Diagonal_MT;
  typedef crhmc_input<MT, NT> INPUT;

  unsigned int _d; // dimension
  MT A;            // matrix A
  VT b;            // vector b, s.t.: Ax=b
  VT lb;
  VT ub;
  int nP;
  // CholObj *solver;
  MT T; // matrix T
  VT y;
  Opts options;
  VT width;
  SimpleBarrier *barrier;
  SpMat Asp; // matrix A
  bool isempty_center = true;
  VT center=VT::Zero(0,1);
  VT w_center;
  int equations() const { return A.rows(); }
  int dimension() const { return A.cols(); }
  int remove_fixed_variables(NT tol = 1e-12) {
    int m = Asp.rows();
    int n = Asp.cols();
    VT d = estimate_width();
    CholObj solver = CholObj(Asp);
    double *w = new double[n];
    for (int i = 0; i < n; i++) {
      w[i] = 1;
    }
    double ac = solver.decompose(w);
    double *out = new double[m];
    double *b_pointer = b.data();
    solver.solve((Tx2 *)b_pointer, (Tx2 *)out);
    VT out_vector =
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(out, m, 1);
    VT x = Asp.transpose() * out_vector;

    x = ((x.array()).abs() < tol).select(0., x);
    std::vector<Triple> freeIndices;
    for (int i = 0; i < n; i++) {
      if (d(i) < tol * (1 + abs(x(i)))) {

      } else {
        freeIndices.push_back(Triple(i, i, 1));
        x(i) = 0.0;
      }
    }
    if (freeIndices.size() != n) {
      SpMat S = SpMat(n, freeIndices.size());
      S.setFromTriplets(freeIndices.begin(), freeIndices.end());

      append_map(S, x);
      return 1;
    }
    return 0;
  }
  int extract_collapsed_variables() {
    SpMat Ac;
    VT bc;
    if (isempty_center) {
      std::tie(center, Ac, bc) = analytic_center(Asp, b, barrier, options);
      isempty_center=false;
    } else {
      std::tie(center, Ac, bc) =
          analytic_center(Asp, b, barrier, options, center);
    }
    if (Ac.rows() == 0) {
      return 0;
    }
    SpMat _A = Asp;
    sparse_stack_v(Ac, _A, Asp);
    b.resize(b.rows() + bc.rows(), 1);
    b << bc, b;

    return 1;
  }
  void updateT() {}
  std::pair<VT, VT> colwiseMinMax(SpMat A) {
    int n = A.cols();
    VT cmax(n);
    VT cmin(n);
    for (int k = 0; k < A.outerSize(); ++k) {
      NT minv = +std::numeric_limits<NT>::infinity();
      NT maxv = -std::numeric_limits<NT>::infinity();
      for (SpMat::InnerIterator it(A, k); it; ++it) {
        minv = std::min(minv, it.value());
        maxv = std::max(maxv, it.value());
      }
      cmin(k) = minv;
      cmax(k) = maxv;
    }
    return std::make_pair(cmin, cmax);
  }
  void nextpow2(VT &a) {
    a = (a.array() == 0).select(1, a);
    a = (((a.array().log()) / std::log(2)).ceil()).matrix();
    a = pow(2, a.array()).matrix();
  }
  std::pair<VT, VT> gmscale(SpMat Asp, NT scltol) {
    int m = Asp.rows();
    int n = Asp.cols();
    SpMat A = Asp.cwiseAbs();
    A.makeCompressed();
    int maxpass = 10;
    NT aratio = 1e+50;
    NT sratio;
    NT damp = 1e-4;
    NT small = 1e-8;
    VT rscale = VT::Ones(m, 1);
    VT cscale = VT::Ones(n, 1);
    VT cmax;
    VT cmin;
    VT rmax;
    VT rmin;
    VT eps = VT::Ones(n, 1) * 1e-12;
    SpMat SA;
    for (int npass = 0; npass < maxpass; npass++) {

      rscale = (rscale.array() == 0).select(1, rscale);
      Diagonal_MT Rinv = (rscale.cwiseInverse()).asDiagonal();
      SA = Rinv * A;
      std::tie(cmin, cmax) = colwiseMinMax(SA);

      //cmin = (cmin + eps).cwiseInverse();
      sratio = (cmax.cwiseQuotient(cmin)).maxCoeff();

      if (npass > 0) {
        cscale = ((cmin.cwiseMax(damp * cmax)).cwiseProduct(cmax)).cwiseSqrt();
      }

      if (npass >= 2 && sratio >= aratio * scltol) {
        break;
      }
      aratio = sratio;
      nextpow2(cscale);
      Diagonal_MT Cinv = (cscale.cwiseInverse()).asDiagonal();
      SA = A * Cinv;
      std::tie(rmin, rmax) = colwiseMinMax(SA.transpose());
      //rmin = (rmin + eps).cwiseInverse();
      rscale = ((rmin.cwiseMax(damp * rmax)).cwiseProduct(rmax)).cwiseSqrt();
      nextpow2(rscale);
    }
    rscale = (rscale.array() == 0).select(1, rscale);
    Diagonal_MT Rinv = (rscale.cwiseInverse()).asDiagonal();
    SA = Rinv * A;
    std::tie(std::ignore, cscale) = colwiseMinMax(SA);
    nextpow2(cscale);
    return std::make_pair(cscale, rscale);
  }

  void rescale(VT x = VT::Zero(0, 1)) {
    if (std::min(equations(), dimension()) <= 1) {
      return;
    }
    VT hess;
    if (x.rows() == 0) {
      hess = VT::Ones(dimension(), 1);
    } else {
      std::tie(std::ignore, hess) = barrier->analytic_center_oracle(x);
      hess = hess + (width.cwiseProduct(width)).cwiseInverse();
    }
    VT scale = (hess.cwiseSqrt()).cwiseInverse();
    SpMat Ain = Asp * scale.asDiagonal();
    VT cscale;
    VT rscale;

    std::tie(cscale, rscale) = gmscale(Ain, 0.9);
    Asp = (rscale.cwiseInverse()).asDiagonal() * Asp;
    b = b.cwiseQuotient(rscale);
    barrier->set_bound(barrier->lb.cwiseProduct(cscale),
                       barrier->ub.cwiseProduct(cscale));
    append_map((cscale.cwiseInverse()).asDiagonal(),VT::Zero(dimension(), 1));
    if (!isempty_center) {
      center = center.cwiseProduct(cscale);
    }
  }
  /*
  std::pair<std::vector<int> , std::vector<int>> nnzPerColumn(){
    std::vector<int> colCounts(n);
    std::vector<int> badCols;

    return std::make_pair(colCounts,badCols)
  }*/


  std::pair<std::vector<int> , std::vector<int>> nnzPerColumn(SpMat A,int threashold){
    int n=A.cols();
    std::vector<int> colCounts(n);
    std::vector<int> badCols;
    for (int k = 0; k < A.outerSize(); ++k) {
      int nnz=0;
        for (SpMat::InnerIterator it(A, k); it; ++it) {
        if(it.value()!=0){
          nnz++;
        }
      }
      colCounts[k]=nnz;
      if(nnz>threashold){
        badCols.push_back(k);
      }
    }
    return std::make_pair(colCounts,badCols);
}
  void splitDenseCols(int maxnz) {
    //std::cout<<"A= \n"<<MT(Asp)<<"\n";
    // Rewrite P so that each cols has no more than maxNZ non-zeros
    int m = Asp.rows();
    int n = Asp.cols();
    if (m <= maxnz) {
      return;
    }
    if (Asp.nonZeros() > maxnz * Asp.cols()) {
      return;
    }
    int numBadCols = 1;

    while (numBadCols > 0) {
      std::vector<int> colCounts(n);
      std::vector<int> badCols;
      numBadCols=0;
      /*
      for (int i = 0; i < Asp.cols(); i++) {
        colCounts[i] = Asp.col(i).nonZeros();
        //std::cout<<i<<" colCounts"<<colCounts[i]<<"\n";
        if (colCounts[i] > maxnz) {
          std::cout<<"Bad Column "<<i<<" colCounts "<<colCounts[i]<<"\n";
          numBadCols++;
          badCols.push_back(i);
        }
      }
      */
      std::tie(colCounts,badCols)=nnzPerColumn(Asp,maxnz);
      if (numBadCols == 0) {
        break;
      }
      m = Asp.rows();
      n = Asp.cols();
      SpMat A_;
      SpMat Aj(m,numBadCols);
      SpMat Ai(numBadCols,n+numBadCols);
      std::vector<Triple> newColumns;
      std::vector<Triple> newRows;
      lb.resize(n + numBadCols);
      ub.resize(n + numBadCols);
      b.resize(m + numBadCols);

      for (int j = 0; j < numBadCols; j++) {
        int i = badCols[j];
        int k = 0;
        for (SpMat::InnerIterator it(Asp, i); it; ++it) {
          if (k >= colCounts[i] / 2) {
        //    std::cout<<"Here "<<it.row()<<" "<<j<<" "<<it.value()<<"\n";
            newColumns.push_back(Triple(it.row(), j, it.value()));
            it.valueRef() = 0;
          }
            k++;
        }
        newRows.push_back(Triple(j, i, 1));
        newRows.push_back(Triple(j, j + n, -1));
        lb(n + j) = lb(i);
        ub(n + j) = ub(i);
        b(m + j) = 0;

      }
      Ai.setFromTriplets(newRows.begin(), newRows.end());
      Aj.setFromTriplets(newColumns.begin(), newColumns.end());
      sparse_stack_h_inplace(Asp, Aj);
      sparse_stack_v(Asp, Ai, A_);
      Asp = A_;
      Asp.makeCompressed();
    }
  }

template <typename MatrixType>
  void append_map(MatrixType S, VT z) {
    b = b - Asp * z;
    Asp = Asp * S;
    y = y + T * z;
    T = T * S;
    updateT();
  }

  void shift_barrier(VT x) {
    int size = x.rows();
    SpMat I = (MT::Identity(size, size)).sparseView();
    append_map(I, x);
    barrier->set_bound(barrier->lb - x, barrier->ub - x);
    if (!isempty_center) {
      center = center - x;
    }
  }
  void reorder() {
    int m = Asp.rows();
    SpMat H;
    H = Asp * SpMat(Asp.transpose()) + MT::Identity(m, m);
    H.makeCompressed();
    Eigen::AMDOrdering<int> ordering;
    PM perm;
    ordering(H, perm);
    H=perm*H*perm.transpose();
    IndexVector m_etree;
    IndexVector firstRowElt;
    Eigen::internal::coletree(H, m_etree, firstRowElt);

    IndexVector post;
    Eigen::internal::treePostorder(int(H.cols()), m_etree, post);

    PM post_perm(m);
    for (int i = 0; i < m; i++)
      post_perm.indices()(i) = post(i);
    perm = perm * post_perm;

    Asp = perm * Asp;
    b = perm * b;
  }
  VT diagL(SpMat &Asp) {
    int m = Asp.cols();
    // SpMat I=(MT::Identity(m, m)).sparseView();
    // SpMat H = SpMat(Asp.transpose()) * Asp + I;
    SpMat H = Asp * SpMat(Asp.transpose());
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                         Eigen::NaturalOrdering<int>>
        cholesky;
    cholesky.analyzePattern(H);
    cholesky.factorize(H);
    SpMat L = cholesky.matrixL();
    return L.diagonal();
  }
  int remove_dependent_rows() {
    int m = Asp.rows();
    VT v = diagL(Asp);
    std::vector<int> indices;
    for (int i = 0; i < m; i++) {
      if ((v(i) > 1e-12) && (v(i) < 1e+64)) {
        indices.push_back(i);
      }
    }
    if (indices.size() == m) {
      return 0;
    }

    Asp = A(indices, Eigen::all).sparseView();
    b = b(indices);
    return 1;
  }
  void simplify() {
    rescale();

    splitDenseCols(options.maxNZ);
    reorder();


    int changed = 1;
    while (changed) {
      while (changed) {
        changed = 0;
        changed += remove_dependent_rows();
        changed += remove_fixed_variables();
        reorder();

      }
      changed += extract_collapsed_variables();
    }
  }

  VT estimate_width() {
    int n = Asp.cols();
    VT hess = VT::Ones(n, 1);
    double *w = new double[n];
    for (int i = 0; i < n; i++) {
      w[i] = 1;
    }
    CholObj solver = CholObj(Asp);
    solver.decompose(w);
    VT w_vector(n, 1);
    solver.leverageScoreComplement((Tx2 *)w_vector.data());
    // VT w_vector = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(w, n,
    w_vector = (w_vector.cwiseMax(0)).cwiseProduct(hess.cwiseInverse());
    VT tau = w_vector.cwiseSqrt();

    delete[] w;
    return tau;
  }
  int dblcmp(NT a, NT b) {
    NT tol = std::numeric_limits<NT>::epsilon();
    return (abs(a - b) < tol * (1 + abs(a) + abs(b)));
  }

  void print() {
    std::cout << "----------------Printing Sparse problem--------------"
              << '\n';
    std::cout << "(m,n) = " << equations() << " , " << dimension() << "\n";
    if(equations() * dimension() > 50){
            std::cout << "too big for complete visulization\n";
            return;
    }
    std::cout << "A=\n";

      std::cout << MT(Asp);
      std::cout << "\n";


    std::cout << "b=\n";
    std::cout << b;
    std::cout << "\n";

    std::cout << "lb=\n";
    std::cout << barrier->lb;
    std::cout << "\n";

    std::cout << "ub=\n";
    std::cout << barrier->ub;
    std::cout << "\n";

    std::cout << "T=\n";
    std::cout << T;
    std::cout << "\n";

    std::cout << "y=\n";
    std::cout << y;
    std::cout << "\n";

    std::cout << "center=\n";
    std::cout << center;
    std::cout << "\n";


  }

crhmcProblem(INPUT &input){
   nP=input.Aeq.cols();
  int nIneq=input.Aineq.rows();
  int nEq=input.Aeq.rows();
  A.resize(nEq+nIneq, nP + nIneq);
  A << input.Aeq, MT::Zero(nEq, nIneq),
        input.Aineq, MT::Identity(nIneq, nIneq);
  b.resize(nEq+nIneq,1);
  b<<input.beq, input.bineq;
  lb.resize(nP+nIneq,1);
  ub.resize(nP+nIneq,1);
  lb<<input.lb, MT::Zero(nIneq, 1);
  ub<<input.ub, MT::Ones(nIneq, 1)* std::numeric_limits<NT>::infinity();

  PreproccessProblem();
}
void PreproccessProblem(){
  int n = dimension();

      /*Move lb=ub to Ax=b*/
      for (int i = 0; i < n; i++) {
        if (dblcmp(lb(i), ub(i))) {
          VT temp = VT::Zero(1, n);
          temp(i) = 1;
          A.conservativeResize(A.rows() + 1, A.cols());
          A.row(A.rows() - 1) = temp;
          b.conservativeResize(b.rows() + 1);
          b(b.rows() - 1) = (lb(i) + ub(i)) / 2;
          lb(i) = -std::numeric_limits<NT>::infinity();
          ub(i) = std::numeric_limits<NT>::infinity();
        }
      }

      barrier = new SimpleBarrier(lb.cwiseMax(-1e7), ub.cwiseMin(1e7));

      Asp = A.sparseView();

      /*Update the transformation Tx + y*/
      T = MT::Zero(nP, n);
      T.block(0, 0, nP, nP) = MT::Identity(nP, nP);
      updateT();
      y = VT::Zero(nP, 1);

      /*Simplify*/
      simplify();

      if (isempty_center) {
        std::tie(center, std::ignore, std::ignore) = analytic_center(Asp, b, barrier, options);
        isempty_center=false;
      }

      shift_barrier(center);
      reorder();

      width = estimate_width();
      if (width.maxCoeff() > 1e9) {
        std::cout << "Domain seems to be unbounded. Either add a Gaussian term "
                     "via f, df, ddf or add bounds to variable via lb and ub."
                  << '\n';
        exit(1);
      }

      //  Recenter again and make sure it is feasible
      VT hess;

      std::tie(center, std::ignore, std::ignore, w_center) =
          lewis_center(Asp, b, barrier, options, center);
      std::tie(std::ignore, hess) = barrier->lewis_center_oracle(center, w_center);
      CholObj solver = CholObj(Asp);
      VT Hinv = hess.cwiseInverse();
      solver.decompose(Hinv.data());
      VT out(equations(), 1);
      VT input = (b - Asp * center);
      solver.solve(input.data(), out.data());
      center = center + (Asp.transpose() * out).cwiseProduct(Hinv);

      if ((center.array() > barrier->ub.array()).any() ||
          (center.array() < barrier->lb.array()).any()) {
        std::cout << "Polytope:Infeasible. The algorithm cannot find a feasible "
                     "point.\n";
        exit(1);
      }

}
  crhmcProblem(HPOLYTOPE &HP) {
    /*Tansform the problem to the form Ax=b lb<=x<=ub*/

   nP = HP.dimension();
    int m = HP.num_of_hyperplanes();
    int n = HP.dimension();

    A.resize(m, n + m);
    A << HP.get_mat(), MT::Identity(m, m);
    b = HP.get_vec();
    n = dimension();
    lb = VT::Zero(n, 1);
    ub = VT::Ones(n) * std::numeric_limits<NT>::infinity();
    PreproccessProblem();
}
};
#endif
