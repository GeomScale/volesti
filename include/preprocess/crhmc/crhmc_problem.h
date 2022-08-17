// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with
// Riemannian Hamiltonian
// Monte Carlo in a Constrained Space"
#ifndef CRHMC_PROBLEM_H
#define CRHMC_PROBLEM_H
#include "Eigen/Eigen"
#include "PackedCSparse/PackedChol.h"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "preprocess/crhmc/analytic_center.h"
#include "preprocess/crhmc/crhmc_input.h"
#include "preprocess/crhmc/crhmc_utils.h"
#include "preprocess/crhmc/lewis_center.h"
#include "preprocess/crhmc/opts.h"
#include "sos/barriers/TwoSidedBarrier.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#ifndef SIMD_LEN
#define SIMD_LEN 0
#endif
const size_t chol_k = (SIMD_LEN == 0) ? 1 : SIMD_LEN;

template <typename Point, typename Input> class crhmc_problem {
public:
  using NT = double;
  using PolytopeType = HPolytope<Point>;
  using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
  using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;
  using IVT = Eigen::Matrix<int, Eigen::Dynamic, 1>;
  using SpMat = Eigen::SparseMatrix<NT>;
  using PM = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;
  using IndexVector = Eigen::Matrix<int, Eigen::Dynamic, 1>;
  using CholObj = PackedChol<chol_k, int>;
  using Triple = Eigen::Triplet<double>;
  using Barrier = TwoSidedBarrier<NT>;
  using Tx = FloatArray<double, chol_k>;
  using Opts = opts<NT>;
  using Diagonal_MT = Eigen::DiagonalMatrix<NT, Eigen::Dynamic>;
  using ArrayXi = Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>;
  using Func = typename Input::Func;
  using Grad = typename Input::Grad;
  using Hess = typename Input::Hess;

  const NT inf = 1e9;

  unsigned int _d; // dimension
  MT A;            // matrix A
  VT b;            // vector b, s.t.: Ax=b
  VT lb;
  VT ub;
  int nP;
  // CholObj *solver;
  SpMat T; // matrix T
  VT y;
  Opts options;
  VT width;
  Barrier barrier;
  SpMat Asp; // matrix A
  bool isempty_center = true;
  VT center = VT::Zero(0, 1);
  VT w_center;
  Hess &ddf;
  Grad &df;
  Func &func;
  bool fZero;     // whether f is completely zero
  bool fHandle;   // whether f is handle or not
  bool dfHandle;  // whether df is handle or not
  bool ddfHandle; // whether ddf is handle or not
  // Assumed each row of T contains at most 1 non-zero
  std::vector<int> Tidx; // T x = x(Tidx) .* Ta
  VT Ta;                 // T x = x(Tidx) .* Ta
  VT Tdf;                // T' * df
  MT Tddf;               // (T.^2)' * ddf
  int equations() const { return Asp.rows(); }
  int dimension() const { return Asp.cols(); }

  VT project(VT const &x) {
    int m = Asp.rows();
    int n = Asp.cols();
    CholObj solver = CholObj(Asp);
    VT hess = barrier.hessian(x);
    VT Hinv = hess.cwiseInverse();
    solver.decompose((Tx *)Hinv.data());
    VT out_vector = VT(m, 1);
    VT in_vector = b.transpose() - Asp * x;
    solver.solve((Tx *)in_vector.data(), (Tx *)out_vector.data());
    out_vector = Asp.transpose() * out_vector;
    return x + (out_vector).cwiseQuotient(hess);
  }

  int remove_fixed_variables(const NT tol = 1e-12) {
    int m = Asp.rows();
    int n = Asp.cols();
    VT d = estimate_width();
    CholObj solver = CholObj(Asp);
    VT w = VT::Ones(n, 1);
    solver.decompose((Tx *)w.data());
    VT out_vector = VT(m, 1);
    solver.solve((Tx *)b.data(), (Tx *)out_vector.data());
    VT x = Asp.transpose() * out_vector;

    x = ((x.array()).abs() < tol).select(0., x);
    std::vector<Triple> freeIndices;
    std::vector<unsigned> indices;
    int nFreeVars = 0;
    for (int i = 0; i < n; i++) {
      if (d(i) < tol * (1 + abs(x(i)))) {
      } else {
        freeIndices.push_back(Triple(i, nFreeVars, 1));
        nFreeVars++;
        indices.push_back(i);
        x(i) = 0.0;
      }
    }

    if (freeIndices.size() != n) {
      SpMat S = SpMat(n, freeIndices.size());
      S.setFromTriplets(freeIndices.begin(), freeIndices.end());
      append_map(S, x);
      barrier.set_bound(barrier.lb(indices), barrier.ub(indices));
      return 1;
    }
    return 0;
  }

  int extract_collapsed_variables() {
    SpMat Ac;
    VT bc;
    if (isempty_center) {
      std::tie(center, Ac, bc) = analytic_center(Asp, b, *this, options);
      isempty_center = false;
    } else {
      std::tie(center, Ac, bc) =
          analytic_center(Asp, b, *this, options, center);
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

  std::pair<VT, VT> colwiseMinMax(SpMat const &A) {
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

  std::pair<VT, VT> gmscale(SpMat &Asp, const NT scltol) {
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

      // cmin = (cmin + eps).cwiseInverse();
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
      // rmin = (rmin + eps).cwiseInverse();
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

  void rescale(const VT x = VT::Zero(0, 1)) {

    if (std::min(equations(), dimension()) <= 1) {
      return;
    }
    VT hess;
    if (x.rows() == 0) {
      hess = VT::Ones(dimension(), 1);
    } else {
      std::tie(std::ignore, hess) = analytic_center_oracle(x);
      hess = hess + (width.cwiseProduct(width)).cwiseInverse();
    }
    VT scale = (hess.cwiseSqrt()).cwiseInverse();
    SpMat Ain = Asp * scale.asDiagonal();
    VT cscale;
    VT rscale;

    std::tie(cscale, rscale) = gmscale(Ain, 0.9);
    Asp = (rscale.cwiseInverse()).asDiagonal() * Asp;
    b = b.cwiseQuotient(rscale);
    barrier.set_bound(barrier.lb.cwiseProduct(cscale),
                      barrier.ub.cwiseProduct(cscale));
    append_map((cscale.cwiseInverse()).asDiagonal(), VT::Zero(dimension(), 1));
    if (!isempty_center) {
      center = center.cwiseProduct(cscale);
    }
  }

  std::pair<std::vector<int>, std::vector<int>>
  nnzPerColumn(SpMat const &A, const int threashold) {
    int n = A.cols();
    std::vector<int> colCounts(n);
    std::vector<int> badCols;
    for (int k = 0; k < A.outerSize(); ++k) {
      int nnz = 0;
      for (SpMat::InnerIterator it(A, k); it; ++it) {
        if (it.value() != 0) {
          nnz++;
        }
      }
      colCounts[k] = nnz;
      if (nnz > threashold) {
        badCols.push_back(k);
      }
    }
    return std::make_pair(colCounts, badCols);
  }

  void splitDenseCols(const int maxnz) {
    //  Rewrite P so that each cols has no more than maxNZ non-zeros
    int m = Asp.rows();
    int n = Asp.cols();
    if (m <= maxnz) {
      return;
    }
    if (Asp.nonZeros() > maxnz * n) {
      return;
    }
    int numBadCols = 1;
    lb = barrier.lb;
    ub = barrier.ub;

    while (numBadCols > 0) {
      m = Asp.rows();
      n = Asp.cols();
      std::vector<int> colCounts(n);
      std::vector<int> badCols;
      numBadCols = 0;
      std::tie(colCounts, badCols) = nnzPerColumn(Asp, maxnz);
      numBadCols = badCols.size();
      if (numBadCols == 0) {
        break;
      }

      SpMat A_;
      SpMat Aj(m, numBadCols);
      SpMat Ai(numBadCols, n + numBadCols);
      std::vector<Triple> newColumns;
      std::vector<Triple> newRows;
      b.conservativeResize(m + numBadCols, 1);
      lb.conservativeResize(n + numBadCols, 1);
      ub.conservativeResize(n + numBadCols, 1);

      for (int j = 0; j < numBadCols; j++) {
        int i = badCols[j];
        int k = 0;
        for (SpMat::InnerIterator it(Asp, i); it; ++it) {
          if (k >= colCounts[i] / 2) {
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
      Asp.prune(0, 0);
      sparse_stack_h_inplace(Asp, Aj);
      sparse_stack_v(Asp, Ai, A_);
      Asp = A_;
    }
    SpMat _T = MT::Zero(T.rows(), ub.rows() - T.cols()).sparseView();
    sparse_stack_h_inplace(T, _T);
    updateT();
    barrier.set_bound(lb, ub);
  }

  template <typename MatrixType>
  void append_map(MatrixType const &S, VT const &z) {
    b = b - Asp * z;
    Asp = Asp * S;
    y = y + T * z;
    T = T * S;
    updateT();
  }

  void shift_barrier(VT const &x) {
    int size = x.rows();
    b = b - Asp * x;
    y = y + T * x;
    barrier.set_bound(barrier.lb - x, barrier.ub - x);
    if (!isempty_center) {
      center = center - x;
    }
  }

  void reorder() {
    if (!options.EnableReordering) {
      return;
    }
    int m = Asp.rows();
    SpMat H;
    H = Asp * SpMat(Asp.transpose()) + MT::Identity(m, m);
    H.makeCompressed();
    Eigen::AMDOrdering<int> ordering;
    PM perm;
    ordering(H, perm);
    // std::cout<<MT(perm)<<"\n\n";
    H = perm * H * perm.transpose();
    IndexVector m_etree;
    IndexVector firstRowElt;
    Eigen::internal::coletree(H, m_etree, firstRowElt);

    IndexVector post;
    Eigen::internal::treePostorder(int(H.cols()), m_etree, post);

    PM post_perm(m);
    for (int i = 0; i < m; i++)
      post_perm.indices()(i) = post(i);
    perm = perm * post_perm;
    // std::cout<<"post_perm= \n"<<MT(perm)<<"\n\n";
    // std::cout<<"m= "<<equations()<<" n= "<<dimension()<<"\n";
    Asp = perm * Asp;
    b = perm * b;
    // std::cout<<"mode \n";
  }

  int remove_dependent_rows(NT tolerance = 1e-12, NT infinity = 1e+64) {
    remove_zero_rows<SpMat, NT>(Asp);
    int m = Asp.rows();
    int n = Asp.cols();
    VT v = VT(m);
    VT w = VT::Ones(n, 1);
    CholObj solver = CholObj(Asp);
    solver.decompose((Tx *)w.data());
    solver.diagL((Tx *)v.data());
    std::vector<int> indices;
    for (int i = 0; i < m; i++) {
      if ((v(i) > tolerance) && (v(i) < infinity)) {
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
#ifdef TIME_KEEPING
    double tstart_rescale, tstart_sparsify, tstart_reorder, tstart_rm_rows,
        tstart_rm_fixed_vars, tstart_ex_colapsed_vars;
    tstart_rescale = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
    rescale();

#ifdef TIME_KEEPING
    std::cout << "Rescale completed in time, ";
    std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_rescale
              << " secs " << std::endl;
#endif
#ifdef TIME_KEEPING
    tstart_sparsify = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
    splitDenseCols(options.maxNZ);
#ifdef TIME_KEEPING
    std::cout << "Split dense columns completed in time, ";
    std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_rescale
              << " secs " << std::endl;
#endif
#ifdef TIME_KEEPING
    tstart_reorder = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
    reorder();
#ifdef TIME_KEEPING
    std::cout << "Reordering completed in time, ";
    std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_reorder
              << " secs " << std::endl;
#endif
    int changed = 1;
    while (changed) {
      while (changed) {
        changed = 0;

#ifdef TIME_KEEPING
        tstart_rm_rows = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
        changed += remove_dependent_rows();
#ifdef TIME_KEEPING
        std::cout << "Removing dependent rows completed in time, ";
        std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_rm_rows
                  << " secs " << std::endl;
#endif

#ifdef TIME_KEEPING
        tstart_rm_fixed_vars = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
        changed += remove_fixed_variables();
#ifdef TIME_KEEPING
        std::cout << "Removing fixed variables completed in time, ";
        std::cout << (double)clock() / (double)CLOCKS_PER_SEC -
                         tstart_rm_fixed_vars
                  << " secs " << std::endl;
#endif
#ifdef TIME_KEEPING
        tstart_reorder = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
        reorder();
#ifdef TIME_KEEPING
        std::cout << "Reordering completed in time, ";
        std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_reorder
                  << " secs " << std::endl;
#endif
      }
#ifdef TIME_KEEPING
      tstart_ex_colapsed_vars = (double)clock() / (double)CLOCKS_PER_SEC;

#endif

      changed += extract_collapsed_variables();
#ifdef TIME_KEEPING
      std::cout << "Extracting collapsed variables completed in time, ";
      std::cout << (double)clock() / (double)CLOCKS_PER_SEC -
                       tstart_ex_colapsed_vars
                << " secs " << std::endl;
#endif
    }
  }

  VT estimate_width() {
    int n = Asp.cols();
    VT hess = VT::Ones(n, 1);
    CholObj solver = CholObj(Asp);
    solver.decompose((Tx *)hess.data());
    VT w_vector(n, 1);
    solver.leverageScoreComplement((Tx *)w_vector.data());
    w_vector = (w_vector.cwiseMax(0)).cwiseProduct(hess.cwiseInverse());
    VT tau = w_vector.cwiseSqrt();

    return tau;
  }

  int doubleVectorEqualityComparison(
      const NT a, const NT b,
      const NT tol = std::numeric_limits<NT>::epsilon()) {
    return (abs(a - b) < tol * (1 + abs(a) + abs(b)));
  }

  void print() {
    std::cerr << "----------------Printing Sparse problem--------------"
              << '\n';
    std::cerr << "(m,n) = " << equations() << " , " << dimension() << "\n";
    if (equations() * dimension() > 50) {
      std::cerr << "too big for complete visulization\n";
      return;
    }
    std::cerr << "A=\n";

    std::cerr << MT(Asp);
    std::cerr << "\n";

    std::cerr << "b=\n";
    std::cerr << b;
    std::cerr << "\n";

    std::cerr << "lb=\n";
    std::cerr << barrier.lb;
    std::cerr << "\n";

    std::cerr << "ub=\n";
    std::cerr << barrier.ub;
    std::cerr << "\n";

    std::cerr << "T=\n";
    std::cerr << MT(T);
    std::cerr << "\n";

    std::cerr << "y=\n";
    std::cerr << y;
    std::cerr << "\n";

    std::cerr << "center=\n";
    std::cerr << center;
    std::cerr << "\n";
  }

  void print(const char *fileName) {
    std::ofstream myfile;
    myfile.open(fileName);
    myfile << Asp.rows() << "  " << Asp.cols() << "\n";

    myfile << MT(Asp);
    myfile << "\n";
    myfile << "\n";

    myfile << b;
    myfile << "\n";
    myfile << "\n";

    myfile << barrier.lb;
    myfile << "\n";
    myfile << "\n";

    myfile << barrier.ub;
    myfile << "\n";
    myfile << "\n";

    myfile << MT(T);
    myfile << "\n";
    myfile << "\n";

    myfile << y;
    myfile << "\n";
    myfile << "\n";

    myfile << center;
  }

  crhmc_problem(Input const &input, Opts _options = Opts())
      : func(input.f), df(input.df), ddf(input.ddf), fHandle(input.fHandle),
        dfHandle(input.dfHandle), ddfHandle(input.ddfHandle),
        fZero(input.fZero) {
    options = _options;
    nP = input.Aeq.cols();
    int nIneq = input.Aineq.rows();
    int nEq = input.Aeq.rows();
    A.resize(nEq + nIneq, nP + nIneq);
    A << input.Aeq, MT::Zero(nEq, nIneq), input.Aineq,
        MT::Identity(nIneq, nIneq);
    b.resize(nEq + nIneq, 1);
    b << input.beq, input.bineq;
    lb.resize(nP + nIneq, 1);
    ub.resize(nP + nIneq, 1);
    lb << input.lb, MT::Zero(nIneq, 1);
    ub << input.ub, MT::Ones(nIneq, 1) * inf;
    Asp.resize(nEq + nIneq, nP + nIneq);
    PreproccessProblem();
  }

  void PreproccessProblem() {
    int n = dimension();

    /*Move lb=ub to Ax=b*/
    for (int i = 0; i < n; i++) {
      if (doubleVectorEqualityComparison(lb(i), ub(i))) {
        VT temp = VT::Zero(1, n);
        temp(i) = 1;
        A.conservativeResize(A.rows() + 1, A.cols());
        A.row(A.rows() - 1) = temp;
        b.conservativeResize(b.rows() + 1);
        b(b.rows() - 1) = (lb(i) + ub(i)) / 2;
        lb(i) = -inf;
        ub(i) = inf;
      }
    }

    barrier.set_bound(lb.cwiseMax(-1e7), ub.cwiseMin(1e7));

    Asp = A.sparseView();
    NT tol = std::numeric_limits<NT>::epsilon();
    Asp.prune(tol, tol);
    /*Update the transformation Tx + y*/
    T = SpMat(nP, n);
    std::vector<Triple> indices;
    for (int i = 0; i < nP; i++) {
      indices.push_back(Triple(i, i, 1));
    }
    T.setFromTriplets(indices.begin(), indices.end());
    Tidx = std::vector<int>(T.rows());
    updateT();
    y = VT::Zero(nP, 1);
    /*Simplify*/
    if (!fZero) {
      fZero = true;
      simplify();
      fZero = false;
    }
    simplify();
#ifdef TIME_KEEPING
    double tstart_rest = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
    if (isempty_center) {
      std::tie(center, std::ignore, std::ignore) =
          analytic_center(Asp, b, *this, options);
      isempty_center = false;
    }
    shift_barrier(center);
#ifdef TIME_KEEPING
    std::cout << "Shift_barrier completed in time, ";
    std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_rest
              << " secs " << std::endl;
#endif
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
#ifdef TIME_KEEPING
    double tstart_find_center = (double)clock() / (double)CLOCKS_PER_SEC;

#endif
    std::tie(center, std::ignore, std::ignore, w_center) =
        lewis_center(Asp, b, *this, options, center);
    std::tie(std::ignore, hess) = lewis_center_oracle(center, w_center);

    CholObj solver = CholObj(Asp);
    VT Hinv = hess.cwiseInverse();
    solver.decompose((Tx *)Hinv.data());
    VT out(equations(), 1);
    VT input = (b - Asp * center);
    solver.solve((Tx *)input.data(), (Tx *)out.data());
    center = center + (Asp.transpose() * out).cwiseProduct(Hinv);
#ifdef TIME_KEEPING
    std::cout << "Finding Center completed in time, ";
    std::cout << (double)clock() / (double)CLOCKS_PER_SEC - tstart_find_center
              << " secs " << std::endl;
#endif
    if ((center.array() > barrier.ub.array()).any() ||
        (center.array() < barrier.lb.array()).any()) {
      std::cout << "Polytope:Infeasible. The algorithm cannot find a feasible "
                   "point.\n";
      exit(1);
    }
  }

  /*Tansform the problem to the form Ax=b lb<=x<=ub*/
  crhmc_problem(PolytopeType const &HP, Opts _options = Opts()) {
    options = _options;
    nP = HP.dimension();
    int m = HP.num_of_hyperplanes();
    int n = HP.dimension();

    A.resize(m, n + m);
    A << HP.get_mat(), MT::Identity(m, m);
    b = HP.get_vec();
    lb.resize(n + m, 1);
    ub.resize(n + m, 1);
    lb = VT(n + m);
    lb << -VT::Ones(n) * inf, VT::Zero(m);
    n = n + m;
    ub = VT::Ones(n) * inf;
    Asp.resize(m, n);
    PreproccessProblem();
  }
  std::pair<VT, VT> analytic_center_oracle(VT const &x) {
    //[~, g, h] = o.f_oracle(x);
    VT g, h;
    std::tie(std::ignore, g, h) = f_oracle(x);
    return std::make_pair(g + barrier.gradient(x), h + barrier.hessian(x));
  }

  std::pair<VT, VT> lewis_center_oracle(VT const &x, VT const &w) {
    //   [~, g, h] = o.f_oracle(x);
    VT g, h;
    std::tie(std::ignore, g, h) = f_oracle(x);
    return std::make_pair(g + w.cwiseProduct(barrier.gradient(x)),
                          h + w.cwiseProduct(barrier.hessian(x)));
  }

  std::tuple<NT, VT, VT> f_oracle(VT x) {
    NT f;
    VT g, h;
    int n = x.rows();
    if (fZero) {
      f = 0;
      g = VT::Zero(n);
      h = VT::Zero(n);
      return std::make_tuple(f, g, h);
    }

    VT z = VT::Zero(n);
    if (fHandle || dfHandle || ddfHandle) {
      z(Tidx, Eigen::all) = Ta.cwiseProduct(x(Tidx, Eigen::all)) + y;
    }
    int k = x.cols();

    if (fHandle) {
      f = func(Point(z));
    } else {
      if (!dfHandle) {
        f = Tdf.transpose() * x;
      } else {
        f = 0;
      }
    }

    if (dfHandle) {
      g = VT::Zero(n, 1);
      g(Tidx, Eigen::all) = Ta.cwiseProduct(df(Point(z)).getCoefficients());
    } else {
      g = Tdf;
    }
    if (ddfHandle) {
      h = VT::Zero(n, 1);
      h(Tidx, Eigen::all) =
          (Ta.cwiseProduct(Ta)).cwiseProduct(ddf(Point(z)).getCoefficients());
    } else {
      h = Tddf;
    }
    return std::make_tuple(f, g, h);
  }

  // By construction each row of T has ar most one nonZero
  void updateT() {
    int n = T.cols();
    int m = T.rows();
    Ta = VT(m);

    for (int k = 0; k < T.outerSize(); ++k) {
      for (SpMat::InnerIterator it(T, k); it; ++it) {
        int pos = (int)it.row();
        int nz = it.col();
        Tidx[pos] = nz;
      }
    }

    Ta = T * VT::Ones(n, 1);
    if (!dfHandle) {
      Tdf = VT::Zero(n);
    }
    if (!ddfHandle) {
      Tddf.resize(n, n);
      Tddf = MT::Zero(n, n);
    }
  }
};
#endif
