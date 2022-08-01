#pragma once
#include <vector>
#include "SparseMatrix.h"

// Problem:
// Compute M = A + B

// Algorithm:
// M = 0
// M(A != 0) += A(A != 0)
// M(B != 0) += B(A != 0)

namespace PackedCSparse {
	template <typename Tx, typename Ti>
	struct AddOutput : SparseMatrix<Tx, Ti>
	{
		UniquePtr<Ti> forwardA;
		UniquePtr<Ti> forwardB;

		template<typename Tx2>
		void initialize(const SparseMatrix<Tx2, Ti>& A, const SparseMatrix<Tx2, Ti>& B)
		{
			pcs_assert(A.initialized() && B.initialized(), "add: bad inputs.");
			pcs_assert(A.n == B.n && A.m == B.m, "add: dimensions mismatch.");

			Ti m = A.m, n = A.n;
			Ti Anz = A.nnz();  Ti* Ap = A.p.get(), * Ai = A.i.get();
			Ti Bnz = B.nnz(); Ti* Bp = B.p.get(), * Bi = B.i.get();
			this->m = A.m; this->n = A.n;

			std::vector<Ti> Ci;
			Ti* Cp = new Ti[n + 1];
			forwardA.reset(new Ti[Anz]);
			forwardB.reset(new Ti[Bnz]);

			Cp[0] = 0;
			for (Ti i = 0; i < n; i++)
			{
				Ti s1 = Ap[i], s2 = Bp[i], end1 = Ap[i + 1], end2 = Bp[i + 1];
				while ((s1 < end1) || (s2 < end2))
				{
					Ti q = Ti(Ci.size());
					Ti i1 = (s1 < end1) ? Ai[s1] : m;
					Ti i2 = (s2 < end2) ? Bi[s2] : m;
					Ti min_i = std::min(i1, i2);
					Ci.push_back(min_i);

					if (i1 == min_i)
						forwardA[s1++] = q;

					if (i2 == min_i)
						forwardB[s2++] = q;
				}
				Cp[i + 1] = Ti(Ci.size());
			}

			this->p.reset(Cp);
			this->i.reset(new Ti[Ci.size()]);
			this->x.reset(pcs_aligned_new<Tx>(Ci.size()));
			std::copy(Ci.begin(), Ci.end(), this->i.get());
		}
	};

	template <typename Tx, typename Ti, typename Tx2 = Tx>
	void add(AddOutput<Tx2, Ti>& o, const SparseMatrix<Tx, Ti>& A, const SparseMatrix<Tx, Ti>& B)
	{
		if (!o.initialized())
			o.initialize(A, B);

		Ti m = o.m, n = o.n;
		Ti Anz = A.nnz(); Ti* Ap = A.p.get(), * Ai = A.i.get(); Tx* Ax = A.x.get();
		Ti Bnz = B.nnz(); Ti* Bp = B.p.get(), * Bi = B.i.get(); Tx* Bx = B.x.get();
		Ti Cnz = o.nnz(); Ti* Cp = o.p.get(), * Ci = o.i.get(); Tx2* Cx = o.x.get();
		Ti* forwardA = o.forwardA.get(), *forwardB = o.forwardB.get();

		for (Ti s = 0; s < Cnz; s++)
			Cx[s] = 0;

		for (Ti s = 0; s < Anz; s++)
			Cx[forwardA[s]] = Ax[s];

		for (Ti s = 0; s < Bnz; s++)
			Cx[forwardB[s]] += Bx[s];
	}

	template <typename Tx, typename Ti, typename Tx2 = Tx>
	AddOutput<Tx2, Ti> add(const SparseMatrix<Tx, Ti>& A, const SparseMatrix<Tx, Ti>& B)
	{
		AddOutput<Tx2, Ti> o;
		add(o, A, B);
		return o;
	}
}