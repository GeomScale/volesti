#pragma once
#include "SparseMatrix.h"

// Problem:
// Compute M = A'

// Algorithm:
// We precompute the mapping from entries of A to entries of At

namespace PackedCSparse {
	template <typename Tx, typename Ti>
	struct TransposeOutput : SparseMatrix<Tx, Ti>
	{
		UniquePtr<Ti> forward;

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx2, Ti>& A)
		{
			pcs_assert(A.initialized(), "transpose: bad inputs.");
			SparseMatrix<Tx, Ti>::initialize(A.n, A.m, A.nnz());

			Ti Am = A.m, An = A.n, * Ap = A.p.get(), * Ai = A.i.get();
			Ti Bm = this->m, Bn = this->n, * Bp = this->p.get(), * Bi = this->i.get();
			Ti nz = A.nnz();

			// compute row counts of A
			Ti* count = new Ti[Bn + 1]();
			
			for (Ti p = 0; p < nz; p++)
				count[Ai[p]]++;

			// compute this->p
			Bp[0] = 0;
			for (Ti i = 0; i < Bn; i++)
			{
				Bp[i + 1] = Bp[i] + count[i];
				count[i] = Bp[i]; // Now, cnt[i] stores the index of the first element in the i-th row
			}

			// compute i and forward
			if (!std::is_same<Tx, bool>::value)
				forward.reset(new Ti[nz]);
			for (Ti j = 0; j < An; j++)
			{
				for (Ti p = Ap[j]; p < Ap[j + 1]; p++)
				{
					Ti q = count[Ai[p]];
					Bi[q] = j;
					if (!std::is_same<Tx, bool>::value)
						forward[p] = q;
					count[Ai[p]]++;
				}
			}

			delete[] count;
		}
	};

	template <typename Tx, typename Ti, typename Tx2 = Tx>
	void transpose(TransposeOutput<Tx2, Ti>& o, const SparseMatrix<Tx, Ti>& A)
	{
		if (!o.initialized())
			o.initialize(A);

		Tx* Ax = A.x.get(); Tx2 *Bx = o.x.get();
		Ti nz = o.nnz(), *forward = o.forward.get();

		if (!std::is_same<Tx2, bool>::value)
		{
			for (Ti s = 0; s < nz; s++)
				Bx[forward[s]] = Tx2(Ax[s]);
		}
	}

	template <typename Tx, typename Ti, typename Tx2 = Tx>
	TransposeOutput<Tx2, Ti> transpose(const SparseMatrix<Tx, Ti>& A)
	{
		TransposeOutput<Tx2, Ti> o;
		transpose(o, A);
		return o;
	}
}
