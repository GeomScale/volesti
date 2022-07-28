#pragma once
#include "SparseMatrix.h"
#include "projinv.h"
#include "outerprod.h"

// Problem:
// Compute M = diag(A' inv(LL') A)

namespace PackedCSparse {
	template<typename Tx, typename Ti>
	struct LeverageOutput : DenseVector<Tx, Ti>
	{
		ProjinvOutput<Tx, Ti> Hinv; // Hinv = inv(H)|_L
		OuterprodOutput<Tx, Ti> tau; // tau = diag(A' Hinv A)

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
		{
			pcs_assert(L.initialized() && A.initialized() && At.initialized(), "leverage: bad inputs.");
			pcs_assert(L.m == L.n && L.n == A.m && L.n == At.n && A.n == At.m, "leverage: dimensions mismatch.");
			DenseVector<Tx, Ti>::initialize(A.n);
			Hinv.initialize(L);
			tau.initialize(A, Hinv, At);
		}
	};

	template<typename Tx, typename Ti, typename Tx2>
	void leverage(LeverageOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
	{
		if (!o.initialized())
			o.initialize(L, A, At);

		Tx T1 = Tx(1.0), T2 = Tx(2.0);
		projinv(o.Hinv, L);

		Ti m = A.m, n = A.n;
		Ti* Sp = o.Hinv.p.get(); Tx* Sv = o.Hinv.x.get();
		for (Ti k = 0; k < m; ++k)
			Sv[Sp[k]] /= T2;

		outerprod(o.tau, A, o.Hinv, At);

		Tx* x = o.x.get(), * tau = o.tau.x.get();
		for (Ti j = 0; j < n; j++)
			x[j] = T2 * tau[j];
	}


	template<typename Tx, typename Ti, typename Tx2>
	LeverageOutput<Tx, Ti> leverage(const SparseMatrix<Tx, Ti>& L, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& At)
	{
		LeverageOutput<Tx, Ti> o;
		leverage(o, L, A, At);
		return o;
	}
}