#pragma once
#include "SparseMatrix.h"

// Problem:
// Compute inv(L L') restricted on L

// Algorithm:
// We need to study this later as this is the bottleneck.
// Document it as a lyx.

namespace PackedCSparse {
	template <typename Tx, typename Ti>
	struct ProjinvOutput : SparseMatrix<Tx, Ti>
	{
		TransposeOutput<bool, Ti> Lt;	// sparsity pattern of the Lt
		UniqueAlignedPtr<Tx> w;			// the row of L we are computing
		UniquePtr<Ti> c;				// c[i] = index the last nonzero on column i in the current L

		void initialize(const SparseMatrix<Tx, Ti>& L)
		{
			pcs_assert(L.initialized(), "chol: bad inputs.");
			pcs_assert(L.n == L.m, "chol: dimensions mismatch.");

			// Copy the sparsity of L
			SparseMatrix<Tx, Ti>::operator=(std::move(L.clone()));

			// allocate workspaces
			Ti n = L.n;
			w.reset(pcs_aligned_new<Tx>(n));
			c.reset(new Ti[n]);
			Lt = transpose<Tx, Ti, bool>(L);
		}
	};

	template <typename Tx, typename Ti>
	void projinv(ProjinvOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& L)
	{
		if (!o.initialized())
			o.initialize(L);

		Tx* Sx = o.x.get(); Ti n = o.n;
		Ti* Li = L.i.get(), * Lp = L.p.get(); Tx* Lv = L.x.get();
		Ti* Lti = o.Lt.i.get(), * Ltp = o.Lt.p.get();
		Tx* w = o.w.get();
		Ti* c = o.c.get();
		Tx T0 = Tx(0), T1 = Tx(1);

		for (Ti k = 0; k < n; k++)
			c[k] = Lp[k + 1] - 1;

		for (Ti k = n - 1; k != -1; k--)
		{
			for (Ti p = Lp[k] + 1; p < Lp[k + 1]; p++)
				w[Li[p]] = Sx[p];

			Tx sum = T1 / Lv[Lp[k]];
			for (Ti p = Ltp[k + 1] - 1; p != Ltp[k] - 1; p--)
			{
				Ti i = Lti[p], Lpi = Lp[i];

				for (Ti q = Lp[i + 1] - 1; q != Lpi; q--)
					fnmadd(sum, Lv[q], w[Li[q]]);
					//sum -= Lv[q] * w[Li[q]];

				sum = sum / Lv[Lpi];
				w[i] = sum;
				Sx[c[i]] = sum;
				c[i]--;
				sum = T0;
			}
		}
	}

	template <typename Tx, typename Ti>
	ProjinvOutput<Tx, Ti> projinv(const SparseMatrix<Tx, Ti>& L)
	{
		ProjinvOutput<Tx, Ti> o;
		projinv(o, L);
		return o;
	}
}