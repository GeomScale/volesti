#pragma once
#include "SparseMatrix.h"

// Problem:
// Compute x = diag(At S Bt)

// Algorithm:
// Note that x = diag(B St A) = grad_H Tr(St A H B)
// We run autodiff on the function Tr(St A H B).
// Hence, the algorithm is essentially same as multiply(A, B) with the same runtime.

namespace PackedCSparse {
	template <typename Tx, typename Ti>
	struct OuterprodOutput : DenseVector<Tx, Ti>
	{
		UniqueAlignedPtr<Tx> s_col;
		UniquePtr<Ti> s_mark;

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx, Ti>& S, const SparseMatrix<Tx2, Ti>& B)
		{
			pcs_assert(A.initialized() && B.initialized() && S.initialized(), "outerprod: bad inputs.");
			pcs_assert(A.m == S.m && S.n == B.n, "outerprod: dimensions mismatch.");

			DenseVector<Tx, Ti>::initialize(A.n);
			s_col.reset(pcs_aligned_new<Tx>(S.m));
			s_mark.reset(new Ti[S.m]);
		}
	};

	template<typename Tx, typename Ti, typename Tx2>
	void outerprod(OuterprodOutput<Tx, Ti>& o, const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx, Ti>& S, const SparseMatrix<Tx2, Ti>& B)
	{
		if (!o.initialized())
			o.initialize(A, S, B);

		Ti Sn = S.n, Sm = S.m, An = A.n;
		Ti* Ap = A.p.get(), * Ai = A.i.get(); Tx2* Ax = A.x.get();
		Ti* Bp = B.p.get(), * Bi = B.i.get(); Tx2* Bx = B.x.get();
		Ti* Sp = S.p.get(), * Si = S.i.get(); Tx* Sx = S.x.get();
		Tx* s_col = o.s_col.get();
		Ti* s_mark = o.s_mark.get();
		Tx* x = o.x.get();

		std::fill(s_mark, s_mark + Sm, Ti(-1));
		std::fill(x, x + An, Tx(0.0));

		for (Ti j = 0; j < Sn; j++)
		{
			for (Ti p = Sp[j]; p < Sp[j + 1]; p++)
			{
				s_col[Si[p]] = Sx[p];
				s_mark[Si[p]] = j;
			}

			for (Ti p = Bp[j]; p < Bp[j + 1]; p++)
			{
				Ti i = Bi[p]; Tx b = Bx[p];
				for (Ti q = Ap[i]; q < Ap[i + 1]; q++)
				{
					Tx a = Ax[q]; Ti a_i = Ai[q];
					if (s_mark[a_i] == j)
					{	//x[i] += s_col[a_i] * a * b;
						fmadd(x[i], s_col[a_i], a * b);
					}
				}
			}
		}
	}

	template <typename Tx, typename Ti, typename Tx2>
	OuterprodOutput<Tx2, Ti> outerprod(const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx, Ti>& S, const SparseMatrix<Tx2, Ti>& B)
	{
		OuterprodOutput<Tx2, Ti> o;
		outerprod(o, A, S, B);
		return o;
	}
}