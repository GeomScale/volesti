#pragma once
#include <memory>
#include <stdexcept>
#include "FloatArray.h"

namespace PackedCSparse {
	static void pcs_assert(bool value, const char* message)
	{
		if (value == false)
			throw std::logic_error(message);
	}

	template <typename T>
	T* pcs_aligned_new(size_t size)
	{
		int alignment = 64; // size of memory cache line
		int offset = alignment - 1 + sizeof(void*);
		void* p1 = (void*)new char[size * sizeof(T) + offset];
		void** p2 = (void**)(((size_t)(p1)+offset) & ~(alignment - 1));
		p2[-1] = p1;
		return (T*)p2;
	}

	template <typename T>
	struct AlignedDeleter
	{
		void operator()(T* p) const
		{
			delete[](char*)(((void**)p)[-1]);
		}
	};

	template<typename T>
	using UniqueAlignedPtr = std::unique_ptr<T[], AlignedDeleter<T>>;

	template<typename T>
	using UniquePtr = std::unique_ptr<T[]>;

	// Tx = Type for entries, Ti = Type for indices.
	// if Tx == bool, the matrix stores only sparsity information
	template <typename Tx, typename Ti>
	struct SparseMatrix
	{
		Ti m = 0;					/* number of rows */
		Ti n = 0;					/* number of columns */
		UniquePtr<Ti> p;			/* column pointers (size n+1) */
		UniquePtr<Ti> i;			/* row indices, size nnz */
		UniqueAlignedPtr<Tx> x;		/* numerical values, size nnz */

		SparseMatrix() = default;

		SparseMatrix(Ti m_, Ti n_, Ti nzmax_)
		{
			initialize(m_, n_, nzmax_);
		}

		bool initialized() const
		{
			return p && i;
		}

		void initialize(Ti m_, Ti n_, Ti nzmax)
		{
			if (nzmax < 1) nzmax = 1;
			m = m_; n = n_;
			p.reset(new Ti[n + 1]);
			i.reset(new Ti[nzmax]);
			if (!std::is_same<Tx, bool>::value)
				x.reset(pcs_aligned_new<Tx>(nzmax));
		}

		Ti nnz() const
		{
			return p[n];
		}

		template <typename Tx2 = Tx, typename Ti2 = Ti>
		SparseMatrix<Tx2, Ti2> clone() const
		{
			SparseMatrix<Tx2, Ti2> C(m, n, nnz());
			Ti* Ap = p.get(), * Ai = i.get(); Tx* Ax = x.get();
			Ti2* Cp = C.p.get(), * Ci = C.i.get(); Tx2* Cx = C.x.get();

			for (Ti s = 0; s <= n; s++)
				Cp[s] = Ti2(Ap[s]);

			Ti nz = nnz();
			for (Ti s = 0; s < nz; s++)
				Ci[s] = Ti2(Ai[s]);

			if (Cx)
			{
				for (Ti s = 0; s < nz; s++)
					Cx[s] = Ax? Tx2(Ax[s]): Tx2(1.0);
			}

			return C;
		}
	};

	template <typename Tx, typename Ti>
	struct DenseVector
	{
		Ti n = 0;					/* number of columns */
		UniqueAlignedPtr<Tx> x;		/* numerical values, size nnz */

		DenseVector() = default;

		DenseVector(Ti n_)
		{
			initialize(n_);
		}

		bool initialized() const
		{
			return bool(x);
		}

		void initialize(Ti n_)
		{
			n = n_;
			x.reset(pcs_aligned_new<Tx>(n_));
		}
	};


	// basic functions
	template <typename Tx, typename Ti>
	SparseMatrix<Tx, Ti> speye(Ti n, Tx* d = nullptr)
	{
		SparseMatrix<Tx, Ti> D(n, n, n);

		for (Ti k = 0; k < n; k++)
		{
			D.i[k] = k;
			D.p[k] = k;
		}
		D.p[n] = n;

		Tx Tx1 = Tx(1.0);
		for (Ti k = 0; k < n; k++)
			D.x[k] = (d ? d[k] : Tx1);
		return D;
	}

	// Solve L out = x
	// Input: L in Tx^{n by n}, x in Tx2^{n}
	// Output: out in Tx2^{n}.
	// If out is provided, we will output to out. Else, output to x.
	template <typename Tx, typename Ti, typename Tx2>
	void lsolve(const SparseMatrix<Tx, Ti>& L, Tx2* x, Tx2* out = nullptr)
	{
		pcs_assert(L.initialized(), "lsolve: bad inputs.");
		pcs_assert(L.n == L.m, "lsolve: dimensions mismatch.");

		Ti n = L.n, * Lp = L.p.get(), * Li = L.i.get(); Tx* Lx = L.x.get();

		if (!out) out = x;
		if (x != out) std::copy(x, x + n, out);

		for (Ti j = 0; j < n; j++)
		{
			Tx2 out_j = out[j] / Lx[Lp[j]];
			out[j] = out_j;

			Ti p_start = Lp[j] + 1, p_end = Lp[j + 1];
			for (Ti p = p_start; p < p_end; p++)
			{	//out[Li[p]] -= Lx[p] * out[j];
				fnmadd(out[Li[p]], out_j, Lx[p]);
			}
		}
	}

	// Solve L' out = x
	// Input: L in Tx^{n by n}, x in Tx2^{n}
	// Output: out in Tx2^{n}.
	// If out is provided, we will output to out. Else, output to x.
	template <typename Tx, typename Ti, typename Tx2>
	void ltsolve(const SparseMatrix<Tx, Ti>& L, Tx2* x, Tx2* out = nullptr)
	{
		pcs_assert(L.initialized(), "ltsolve: bad inputs.");
		pcs_assert(L.n == L.m, "ltsolve: dimensions mismatch.");

		Ti n = L.n, * Lp = L.p.get(), * Li = L.i.get(); Tx* Lx = L.x.get();

		if (!out) out = x;
		if (x != out) std::copy(x, x + n, out);

		for (Ti j = n - 1; j != -1; j--)
		{
			Tx2 out_j = out[j];

			Ti p_start = Lp[j] + 1, p_end = Lp[j + 1];
			for (Ti p = p_start; p < p_end; p++)
			{   //out[j] -= Lx[p] * out[Li[p]];
				fnmadd(out_j, out[Li[p]], Lx[p]);
			}

			out[j] = out_j / Tx2(Lx[Lp[j]]);
		}
	}

	// Update y <-- y + A x
	// Input: A in Tx^{n by n}, x, y in Tx2^{n}
	template <typename Tx, typename Ti, typename Tx2>
	void gaxpy(const SparseMatrix<Tx, Ti>& A, const Tx2* x, Tx2* y)
	{
		pcs_assert(A.initialized(), "gaxpy: bad inputs.");
		Ti m = A.m, n = A.n, * Ap = A.p.get(), * Ai = A.i.get(); Tx* Ax = A.x.get();

		for (Ti j = 0; j < n; j++)
		{
			Tx2 x_j = x[j];

			Ti p_start = Ap[j], p_end = Ap[j + 1];
			for (Ti p = p_start; p < p_end; p++)
			{   //y[Ai[p]] += Ax[p] * x[j];
				fmadd(y[Ai[p]], x_j, Ax[p]);
			}
		}
	}
};
