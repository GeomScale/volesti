#pragma once
#include <vector>
#include <algorithm>
#include "SparseMatrix.h"

// Problem:
// Compute M = A diag(w) B

// Algorithm:
// Compute M col by col

namespace PackedCSparse {
	template <typename Tx, typename Ti>
	struct MultiplyOutput : SparseMatrix<Tx, Ti>
	{
		UniqueAlignedPtr<Tx> c;

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& B)
		{
			pcs_assert(A.initialized() && B.initialized(), "multiply: bad inputs.");
			pcs_assert(A.n == B.m, "multiply: dimensions mismatch.");

			Ti m = A.m, n = B.n;
			Ti* Ap = A.p.get(), * Ai = A.i.get();
			Ti* Bp = B.p.get(), * Bi = B.i.get();

			this->c.reset(pcs_aligned_new<Tx>(m));

			Ti* last_j = new Ti[m];
			for (Ti i = 0; i < m; i++)
			{
				last_j[i] = -1;
				this->c[i] = Tx(0.0);
			}

			Ti* Cp = new Ti[size_t(n)+1];
			std::vector<Ti> Ci;

			Cp[0] = 0;
			for (Ti j1 = 0; j1 < n; j1++)
			{
				for (Ti p1 = Bp[j1]; p1 < Bp[j1 + 1]; p1++)
				{
					Ti j2 = Bi[p1];
					for (Ti p2 = Ap[j2]; p2 < Ap[j2 + 1]; p2++)
					{
						Ti i = Ai[p2];
						if (last_j[i] != j1)
						{
							last_j[i] = j1;
							Ci.push_back(i);
						}
					}
				}
				Cp[j1 + 1] = Ti(Ci.size());
			}
			delete[] last_j;

			for (Ti j = 0; j < n; j++)
				std::sort(Ci.begin() + Cp[j], Ci.begin() + Cp[j + 1]);

			this->m = m; this->n = n;
			this->x.reset(pcs_aligned_new<Tx>(Ci.size()));
			this->p.reset(Cp);
			this->i.reset(new Ti[Ci.size()]);
			std::copy(Ci.begin(), Ci.end(), this->i.get());
		}
	};

	template <typename Tx, typename Ti, typename Tx2, bool has_weight>
	void multiply_general(MultiplyOutput<Tx2, Ti>& o, const SparseMatrix<Tx, Ti>& A, const Tx2* w, const SparseMatrix<Tx, Ti>& B)
	{
		if (!o.initialized())
			o.initialize(A, B);

		Ti m = o.m, n = o.n;
		Ti* Ap = A.p.get(), * Ai = A.i.get(); Tx* Ax = A.x.get();
		Ti* Bp = B.p.get(), * Bi = B.i.get(); Tx* Bx = B.x.get();
		Ti* Cp = o.p.get(), * Ci = o.i.get(); Tx2* Cx = o.x.get();
		Tx2* c = o.c.get(); // initialized to 0

		const Tx2 T0 = Tx2(0);
		for (Ti j1 = 0; j1 < n; j1++)
		{
			for (Ti p1 = Bp[j1]; p1 < Bp[j1 + 1]; p1++)
			{
				Ti j2 = Bi[p1];
				Tx2 beta = has_weight? (Tx2(Bx[p1]) * w[j2]) : Tx2(Bx[p1]);

				for (Ti p2 = Ap[j2]; p2 < Ap[j2 + 1]; p2++)
				{
					//x[Ai[p2]] += beta * Ax[p2];
					fmadd(c[Ai[p2]], beta, Ax[p2]);
				}
			}

			for (Ti p1 = Cp[j1]; p1 < Cp[j1 + 1]; p1++)
			{
				Cx[p1] = c[Ci[p1]];
				c[Ci[p1]] = T0; // ensure c is 0 after the call
			}
		}
	}

	template <typename Tx, typename Ti>
	void multiply(MultiplyOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& A, const SparseMatrix<Tx, Ti>& B)
	{
		multiply_general<Tx, Ti, Tx, false>(o, A, nullptr, B);
	}

	template <typename Tx, typename Ti, typename Tx2>
	void multiply(MultiplyOutput<Tx2, Ti>& o, const SparseMatrix<Tx, Ti>& A, const Tx2* w, const SparseMatrix<Tx, Ti>& B)
	{
		multiply_general<Tx, Ti, Tx2, true>(o, A, w, B);
	}

	template <typename Tx, typename Ti, typename Tx2>
	MultiplyOutput<Tx2, Ti> multiply(const SparseMatrix<Tx, Ti>& A, const Tx2* w, const SparseMatrix<Tx, Ti>& B)
	{
		MultiplyOutput<Tx2, Ti> o;
		multiply(o, A, w, B);
		return o;
	}

	template <typename Tx, typename Ti>
	MultiplyOutput<Tx, Ti> multiply(const SparseMatrix<Tx, Ti>& A, const SparseMatrix<Tx, Ti>& B)
	{
		MultiplyOutput<Tx, Ti> o;
		multiply(o, A, B);
		return o;
	}
}
