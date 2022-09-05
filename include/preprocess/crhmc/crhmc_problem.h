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
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#ifndef SIMD_LEN
#define SIMD_LEN 0
#endif
const size_t chol_k = (SIMD_LEN == 0) ? 1 : SIMD_LEN;
template <typename Point, typename Input>
class crhmc_problem {
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
  using Triple = Eigen::Triplet<NT>;
  using Barrier = TwoSidedBarrier<Point>;
  using Tx = FloatArray<double, chol_k>;
  using Opts = opts<NT>;
  using Diagonal_MT = Eigen::DiagonalMatrix<NT, Eigen::Dynamic>;
  using Func = typename Input::Func;
  using Grad = typename Input::Grad;
  using Hess = typename Input::Hess;

  unsigned int _d; // dimension
  // Problem variables Ax=b st lb<=x<=ub
  MT A;            // matrix A input matrix
  SpMat Asp;       // problem matrix A in Sparse form
  VT b;            // vector b, s.t.: Ax=b
  VT lb;           // Lower bound for output coordinates
  VT ub;           // Upper bound for output coordinates
  Barrier barrier; // Class that holds functions that handle the log barrier for
                   // lb and ub
  Opts options;    // problem parameters
  // Transformation (T,y) such that the new variables x
  // can be tranform to the original z (z= Tx+y)
  SpMat T;
  VT y;
  // Non zero indices and values for fast tranform
  std::vector<int> Tidx; // T x = x(Tidx) .* Ta
  VT Ta;                 // T x = x(Tidx) .* Ta
  bool isempty_center = true;
  VT center = VT::Zero(0, 1); // Resulting polytope Lewis or Analytic center
  VT w_center;// weights of the lewis center

  VT width; // width of the varibles
  int nP;//input dimension

  Func &func;     // function handle
  Grad &df;       // gradient handle
  Hess &ddf;      // hessian handle
  bool fZero;     // whether f is completely zero
  bool fHandle;   // whether f is handle or not
  bool dfHandle;  // whether df is handle or not
  bool ddfHandle; // whether ddf is handle or not
#ifdef TIME_KEEPING
//Timing information
  std::chrono::duration<double> rescale_duration, sparsify_duration,
      reordering_duration, rm_rows_duration, rm_fixed_vars_duration,
      ex_collapsed_vars_duration, shift_barrier_duration, lewis_center_duration;
#endif
  const NT inf = options.max_coord + 1; // helper for barrier handling
  int equations() const { return Asp.rows(); }
  int dimension() const { return Asp.cols(); }

  // Remove varibles that have width under some tolerance
  int remove_fixed_variables(const NT tol = 1e-12) {
    int m = Asp.rows();
    int n = Asp.cols();
    VT d = estimate_width();
    CholObj solver = CholObj(transform_format<SpMat,NT,int>(Asp));
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
  // Rescale the polytope for numerical stability
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

    std::tie(cscale, rscale) = gmscale<SpMat, VT, NT>(Ain, 0.9);
    Asp = (rscale.cwiseInverse()).asDiagonal() * Asp;
    b = b.cwiseQuotient(rscale);
    barrier.set_bound(barrier.lb.cwiseProduct(cscale),
                      barrier.ub.cwiseProduct(cscale));
    append_map((cscale.cwiseInverse()).asDiagonal(), VT::Zero(dimension(), 1));
    if (!isempty_center) {
      center = center.cwiseProduct(cscale);
    }
  }

  //  Rewrite P so that each cols has no more than maxNZ non-zeros
  void splitDenseCols(const int maxnz) {
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
    //until there are columns with more than maxNZ elements
    while (numBadCols > 0) {
      m = Asp.rows();
      n = Asp.cols();
      std::vector<int> colCounts(n);
      std::vector<int> badCols;
      numBadCols = 0;
      //find the columns with count larger than maxNZ
      std::tie(colCounts, badCols) = nnzPerColumn(Asp, maxnz);
      numBadCols = badCols.size();
      if (numBadCols == 0) {
        break;
      }
      //create a new variable for each one and update Asp, b, lb, ub, T, y accordingly
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
  // Change A and the correpsonding Transformation
  template <typename MatrixType>
  void append_map(MatrixType const &S, VT const &z) {
    b = b - Asp * z;
    Asp = Asp * S;
    y = y + T * z;
    T = T * S;
    updateT();
  }
  // Shift the problem with respect to x
  void shift_barrier(VT const &x) {
    int size = x.rows();
    b = b - Asp * x;
    y = y + T * x;
    barrier.set_bound(barrier.lb - x, barrier.ub - x);
    if (!isempty_center) {
      center = center - x;
    }
  }
  // Reorder the polytope accordint to the AMD Reordering for better sparsity
  // pattern in the Cholesky decomposition
  void reorder() {
    if (!options.EnableReordering) {
      return;
    }
    Asp.prune(0.0);
    Asp.makeCompressed();
    int m = Asp.rows();
    SpMat H;
    H = Asp * SpMat(Asp.transpose()) + MT::Identity(m, m);
    H.makeCompressed();
    PM permed = permuteMatAMD(H);
    H = permed * H * permed.transpose();
    PM post_perm = postOrderPerm(H);
    PM perm = permed * post_perm;
    Asp = perm * Asp;
    b = perm * b;
  }
//Using the Cholesky decomposition remove dependent rows in the systen Asp*x=b
  int remove_dependent_rows(NT tolerance = 1e-12, NT infinity = 1e+64) {
    //this approach does not work with 0 collumns
    remove_zero_rows<SpMat, NT>(Asp);
    int m = Asp.rows();
    int n = Asp.cols();
    VT v = VT(m);
    VT w = VT::Ones(n, 1);
    CholObj solver = CholObj(transform_format<SpMat,NT,int>(Asp));
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
//Apply number of operations that simplify the problem
  void simplify() {
#ifdef TIME_KEEPING
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::system_clock::now();
#endif
    rescale();
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    rescale_duration += end - start;
    start = std::chrono::system_clock::now();
#endif
    splitDenseCols(options.maxNZ);
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    sparsify_duration += end - start;
    start = std::chrono::system_clock::now();
#endif
    reorder();
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    reordering_duration += end - start;
#endif
    int changed = 1;
    while (changed) {
      while (changed) {
        changed = 0;
#ifdef TIME_KEEPING
        start = std::chrono::system_clock::now();
#endif
        changed += remove_dependent_rows();
#ifdef TIME_KEEPING
        end = std::chrono::system_clock::now();
        rm_rows_duration += end - start;
        start = std::chrono::system_clock::now();
#endif
        changed += remove_fixed_variables();
#ifdef TIME_KEEPING
        end = std::chrono::system_clock::now();
        rm_fixed_vars_duration += end - start;
        start = std::chrono::system_clock::now();
#endif
        reorder();
#ifdef TIME_KEEPING
        end = std::chrono::system_clock::now();
        reordering_duration += end - start;
#endif
      }
#ifdef TIME_KEEPING
      start = std::chrono::system_clock::now();
#endif
      changed += extract_collapsed_variables();
#ifdef TIME_KEEPING
      end = std::chrono::system_clock::now();
      ex_collapsed_vars_duration += end - start;
#endif
    }
  }

  VT estimate_width() {
    int n = Asp.cols();
    VT hess = VT::Ones(n, 1);
    CholObj solver = CholObj(transform_format<SpMat,NT,int>(Asp));
    solver.decompose((Tx *)hess.data());
    VT w_vector(n, 1);
    solver.leverageScoreComplement((Tx *)w_vector.data());
    w_vector = (w_vector.cwiseMax(0)).cwiseProduct(hess.cwiseInverse());
    VT tau = w_vector.cwiseSqrt();

    return tau;
  }

  void print(std::string const message = "Printing Sparse problem") {
    std::cerr << "----------------" << message << "--------------" << '\n';
    std::cerr << "(m,n) = " << equations() << " , " << dimension()
              << " nnz= " << Asp.nonZeros() << "\n";
    if (equations() > 20 || dimension() > 20) {
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
//Class constructor
  crhmc_problem(Input const &input, Opts _options = Opts())
      : options(_options), func(input.f), df(input.df), ddf(input.ddf),
        fZero(input.fZero), fHandle(input.fHandle), dfHandle(input.dfHandle),
        ddfHandle(input.ddfHandle) {
#ifdef TIME_KEEPING
    rescale_duration = sparsify_duration = reordering_duration =
        rm_rows_duration = rm_fixed_vars_duration = ex_collapsed_vars_duration =
            shift_barrier_duration = lewis_center_duration =
                std::chrono::duration<double>::zero();
#endif

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
  // Initialization funciton
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
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::system_clock::now();
#endif
    if (isempty_center) {
      std::tie(center, std::ignore, std::ignore) =
          analytic_center(Asp, b, *this, options);
      isempty_center = false;
    }
    shift_barrier(center);
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    shift_barrier_duration += end - start;
#endif
    reorder();

    width = estimate_width();
    if (width.maxCoeff() > 1e9) {
      std::cerr << "Domain seems to be unbounded. Either add a Gaussian term "
                   "via f, df, ddf or add bounds to variable via lb and ub."
                << '\n';
      exit(1);
    }
    //  Recenter again and make sure it is feasible
    VT hess;
#ifdef TIME_KEEPING
    start = std::chrono::system_clock::now();
#endif
    std::tie(center, std::ignore, std::ignore, w_center) =
        lewis_center(Asp, b, *this, options, center);
    std::tie(std::ignore, hess) = lewis_center_oracle(center, w_center);
    CholObj solver = CholObj(transform_format<SpMat,NT,int>(Asp));
    VT Hinv = hess.cwiseInverse();
    solver.decompose((Tx *)Hinv.data());
    VT out(equations(), 1);
    VT input = (b - Asp * center);
    solver.solve((Tx *)input.data(), (Tx *)out.data());
    center = center + (Asp.transpose() * out).cwiseProduct(Hinv);
#ifdef TIME_KEEPING
    end = std::chrono::system_clock::now();
    lewis_center_duration += end - start;
#endif
    if ((center.array() > barrier.ub.array()).any() ||
        (center.array() < barrier.lb.array()).any()) {
      std::cerr << "Polytope:Infeasible. The algorithm cannot find a feasible "
                   "point.\n";
      exit(1);
    }
#ifdef TIME_KEEPING
    std::cerr << "Rescale completed in time, ";
    std::cerr << rescale_duration.count() << " secs " << std::endl;
    std::cerr << "Split dense columns completed in time, ";
    std::cerr << sparsify_duration.count() << " secs " << std::endl;
    std::cerr << "Reordering completed in time, ";
    std::cerr << reordering_duration.count() << " secs " << std::endl;
    std::cerr << "Removing dependent rows completed in time, ";
    std::cerr << rm_rows_duration.count() << " secs " << std::endl;
    std::cerr << "Removing fixed variables completed in time, ";
    std::cerr << rm_fixed_vars_duration.count() << " secs " << std::endl;
    std::cerr << "Extracting collapsed variables completed in time, ";
    std::cerr << ex_collapsed_vars_duration.count() << " secs " << std::endl;
    std::cerr << "Shift_barrier completed in time, ";
    std::cerr << shift_barrier_duration.count() << " secs " << std::endl;
    std::cerr << "Finding Center completed in time, ";
    std::cerr << lewis_center_duration.count() << " secs " << std::endl;
#endif
  }

  // Gradient and hessian of for the analytic center
  std::pair<VT, VT> analytic_center_oracle(VT const &x) {
    VT g, h;
    std::tie(std::ignore, g, h) = f_oracle(x);
    return std::make_pair(g + barrier.gradient(x), h + barrier.hessian(x));
  }
  // Gradient and hessian of for the lewis center
  std::pair<VT, VT> lewis_center_oracle(VT const &x, VT const &w) {
    VT g, h;
    std::tie(std::ignore, g, h) = f_oracle(x);
    return std::make_pair(g + w.cwiseProduct(barrier.gradient(x)),
                          h + w.cwiseProduct(barrier.hessian(x)));
  }
  // Function that uses the transformation (T,y) to apply the function to the
  // original variables
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
    // Take the correpsonding point in the original space
    VT z = VT::Zero(n);
    if (fHandle || dfHandle || ddfHandle) {
      z(Tidx, Eigen::all) = Ta.cwiseProduct(x(Tidx, Eigen::all)) + y;
    }

    // If the function is given evaluate it at the original point
    if (fHandle) {
      f = func(Point(z));
    } else {
      f = 0;
    }
    // If the gradient is given evaluate it at the original point
    if (dfHandle) {
      g = VT::Zero(n, 1);
      g(Tidx, Eigen::all) = Ta.cwiseProduct(df(Point(z)).getCoefficients());
    } else {
      g = VT::Zero(n, 1);
    }
    // If the hessian is given evaluate it at the original point
    if (ddfHandle) {
      h = VT::Zero(n, 1);
      h(Tidx, Eigen::all) =
          (Ta.cwiseProduct(Ta)).cwiseProduct(ddf(Point(z)).getCoefficients());
    } else {
      h = VT::Zero(n, 1);
    }
    return std::make_tuple(f, -g, h);
  }

  // Update the indices and values vectors of the matrix T
  void updateT() {
    int n = T.cols();
    int m = T.rows();
    Ta = VT(m);
    // By construction each row of T has ar most one nonZero
    for (int k = 0; k < T.outerSize(); ++k) {
      for (SpMat::InnerIterator it(T, k); it; ++it) {
        int pos = (int)it.row();
        int nz = it.col();
        Tidx[pos] = nz;
      }
    }

    Ta = T * VT::Ones(n, 1);
  }
};
#endif
