// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis
// Copyright (c) 2022 Ioannis Iakovidis

// This file is converted from PolytopeSamplerMatlab
//(https://github.com/ConstrainedSampler/PolytopeSamplerMatlab/blob/master/code/solver/PackedCSparse/PackedChol.h) by Ioannis Iakovidis

#pragma once
#include "Eigen/Eigen"
#include "SparseMatrix.h"
#include "chol.h"
#include "leverage.h"
#include "leverageJL.h"
#include "multiply.h"
#include "qd/dd_real.h"
#include <random>
#include <vector>
using namespace PackedCSparse;

template <typename Tin, typename Tout>
void get_slice(Tout *out, Tin *in, size_t n, size_t idx) {
  for (size_t j = 0; j < n; j++)
    out[j] = to_double(get(in[j], idx));
}

template <typename Tin, typename Tout>
void set_slice(Tout *out, Tin *in, size_t n, size_t idx) {
  for (size_t j = 0; j < n; j++)
    set(out[j], idx, to_double(in[j]));
}

template <int k, typename Ti> struct PackedChol {
  using Tx = double;
  using Tx2 = FloatArray<double, k>;
  using Te = dd_real;

  // parameters
  SparseMatrix<Tx, Ti> A;
  SparseMatrix<Tx, Ti> At;
  UniqueAlignedPtr<Tx2> w;
  Tx accuracyThreshold = 1e-6;
  std::vector<size_t>
      exactIdx; // k size array. Indices we perform high precision calculation
  std::vector<size_t>
      numExact; // number of times we perform high precision decompose (length
                // k+1, the last one records how many times we do decompose)
  bool decomposed = false;

  // preprocess info for different CSparse operations (PackedDouble)
  MultiplyOutput<Tx2, Ti> H;         // cache for H = A W A'
  CholOutput<Tx2, Ti> L;             // cache for L = chol(H)
  LeverageOutput<Tx2, Ti> diagP;     // cache for L = chol(H)
  LeverageJLOutput<Tx2, Ti> diagPJL; // cache for L = chol(H)

  // preprocess info for different CSparse operations (dd_real)
  MultiplyOutput<Te, Ti> H_exact;         // cache for H = A W A'
  CholOutput<Te, Ti> L_exact;             // cache for L = chol(H)
  LeverageOutput<Te, Ti> diagP_exact;     // cache for L = chol(H)
  LeverageJLOutput<Te, Ti> diagPJL_exact; // cache for L = chol(H)
  SparseMatrix<Te, Ti> Le[k];             // store output of L_exact

  PackedChol(const SparseMatrix<Tx, Ti> &A_) {
    A = std::move(A_.clone());
    At = transpose(A);
    w.reset(pcs_aligned_new<Tx2>(A.n));
    numExact.resize(k + 1);
  }

  void setSeed(unsigned long long seed) {
    diagPJL.gen.seed(seed);
    diagPJL_exact.gen.seed(seed);
  }

  bool allExact() { return exactIdx.size() == k; }

  bool hasExact() { return exactIdx.size() > 0; }

  template <typename Tv2_> Tx2 decompose(const Tv2_ *w_in) {

    Tx2 acc = Tx2(0.0);

    // record w
    Ti n = A.n;

    for (Ti j = 0; j < n; j++) {
      w[j] = w_in[j];
    }
    // compute chol
    ++numExact[k];
    if (accuracyThreshold > 0.0 ||
        !decomposed) // the first time we call, always run the double chol.
    {
      multiply(H, A, w.get(), At);
      chol(L, H);
      decomposed = true;

      exactIdx.clear();
      acc = estimateAccuracy();
      for (size_t i = 0; i < k; i++) {
        if (get(acc, i) >=
            accuracyThreshold) // >= is important for the case accuracyThreshold
                               // = 0.0, we need to compute everything exactly
          exactIdx.push_back(i);
      }
    } else if (!allExact()) {
      exactIdx.clear();
      for (size_t i = 0; i < k; i++)
        exactIdx.push_back(i);
    }

    if (hasExact()) {
      Te *w_exact = new Te[n];

      for (size_t i : exactIdx) {
        ++numExact[i];
        get_slice(w_exact, w.get(), n, i);
        multiply(H_exact, A, w_exact, At);
        chol(L_exact, H_exact);

        // copy result to Le[i]
        if (!Le[i].initialized())
          Le[i] = std::move(L_exact.template clone<Te, Ti>());
        else {
          Ti nz = L_exact.nnz();
          for (Ti s = 0; s < nz; ++s)
            Le[i].x[s] = (L_exact.x[s]);
        }
      }

      delete[] w_exact;
    }
    return acc;
  }

  Tx2 logdet() {
    pcs_assert(decomposed, "logdet: Need to call decompose first.");

    Ti m = A.m;
    Tx2 ret = Tx2(0);

    if (!allExact()) {
      Ti *Lp = L.p.get();
      Tx2 *Lx = L.x.get();
      for (Ti j = 0; j < m; j++)
        ret += log(Lx[Lp[j]]);
    }

    if (hasExact()) {
      for (size_t i : exactIdx) {
        Te ret_e = 0.0;
        Ti *Lp = Le[i].p.get();
        Te *Lx = Le[i].x.get();

        for (Ti j = 0; j < m; j++)
          ret_e += log(Lx[Lp[j]]);

        set(ret, i, to_double(ret_e));
      }
    }

    return ret * Tx2(2.0);
  }

  void diagL(Tx2 *out) {
    pcs_assert(decomposed, "diagL: Need to call decompose first.");

    Ti m = A.m;

    if (!allExact()) {
      Ti *Li = L.i.get(), *Lp = L.p.get();
      Tx2 *Lx = L.x.get();
      for (Ti j = 0; j < m; j++)
        out[j] = Lx[Lp[j]];
    }

    if (hasExact()) {
      for (size_t i : exactIdx) {
        Ti *Lp = Le[i].p.get();
        Te *Lx = Le[i].x.get();

        for (Ti j = 0; j < m; j++)
          set(out[j], i, to_double(Lx[Lp[j]]));
      }
    }
  }

  SparseMatrix<double, Ti> getL(Ti i) {
    pcs_assert(decomposed, "getL: Need to call decompose first.");

    Ti m = L.m, n = L.n, nz = L.nnz();
    SparseMatrix<double, Ti> out(m, n, nz);

    Ti *outp = out.p.get(), *Lp = L.p.get();
    Ti *outi = out.i.get(), *Li = L.i.get();

    for (Ti s = 0; s <= n; s++)
      outp[s] = Lp[s];

    for (Ti s = 0; s < nz; s++)
      outi[s] = Li[s];

    bool isExact = false;
    for (size_t i_ : exactIdx) {
      if (i_ == i)
        isExact = true;
    }

    double *outx = out.x.get();
    if (isExact) {
      Te *Lx = Le[i].x.get();
      for (Ti s = 0; s < nz; s++)
        outx[s] = to_double(Lx[s]);
    } else {
      Tx2 *Lx = L.x.get();
      for (Ti s = 0; s < nz; s++)
        outx[s] = get(Lx[s], i);
    }

    return std::move(out);
  }

  void solve(Tx2 *b, Tx2 *out) {
    pcs_assert(decomposed, "solve: Need to call decompose first.");

    if (!allExact()) {
      lsolve(L, b, out);
      ltsolve(L, out, out);
    }

    if (hasExact()) {
      Ti m = A.m;
      Te *b_exact = new Te[m];
      Te *out_exact = new Te[m];

      for (size_t i : exactIdx) {
        get_slice(b_exact, b, m, i);
        lsolve(Le[i], b_exact, out_exact);
        ltsolve(Le[i], out_exact, out_exact);
        set_slice(out, out_exact, m, i);
      }

      delete[] b_exact;
      delete[] out_exact;
    }
  };

  void leverageScoreComplement(Tx2 *out) {
    pcs_assert(decomposed,
               "leverageScoreComplement: Need to call decompose first.");

    Ti n = A.n, m = A.m;

    if (!allExact()) {
      Tx2 T1 = Tx2(1.0), T2 = Tx2(2.0);
      leverage(diagP, L, A, At);

      Tx2 *tau = diagP.x.get();
      for (Ti j = 0; j < n; j++)
        out[j] = T1 - tau[j] * w[j];
    }

    if (hasExact()) {
      Te T1 = Te(1.0), T2 = Te(2.0);
      for (size_t i : exactIdx) {
        leverage(diagP_exact, Le[i], A, At);

        Te *tau = diagP_exact.x.get();
        for (Ti j = 0; j < n; j++)
          set(out[j], i, to_double(T1 - tau[j] * get(w[j], i)));
      }
    }
  }

  void leverageScoreComplementJL(Tx2 *out, size_t JL_k) {
    pcs_assert(decomposed,
               "leverageScoreComplementJL: Need to call decompose first.");

    Ti m = A.m, n = A.n;

    if (!allExact()) {
      Tx2 T1 = Tx2(1.0), T2 = Tx2(2.0);
      leverageJL(diagPJL, L, A, At, JL_k);

      Tx2 *tau = diagPJL.x.get();
      for (Ti j = 0; j < n; j++)
        out[j] = T1 - tau[j] * w[j];
    }

    if (hasExact()) {
      Te T1 = Te(1.0), T2 = Te(2.0);
      for (size_t i : exactIdx) {
        leverageJL(diagPJL_exact, Le[i], A, At, JL_k);

        Te *tau = diagPJL_exact.x.get();
        for (Ti j = 0; j < n; j++)
          set(out[j], i, to_double(T1 - tau[j] * get(w[j], i)));
      }
    }
  }

  Tx2 estimateAccuracy() {
    pcs_assert(decomposed, "estimateAccuracy: Need to call decompose first.");

    return cholAccuracy(diagPJL, L, A, At, w.get());
  }
};
