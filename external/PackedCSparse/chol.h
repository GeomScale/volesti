#pragma once
#include <vector>
#include <queue>
#include "SparseMatrix.h"
#include "transpose.h"

// Problem:
// Compute chol(A)

// Algorithm:
// We need to study this later as this is the bottleneck.
// Document it as a lyx.
// chol_up_looking:
//		Compute L row by row
//		This is faster when it is compute bound.
//
// chol_left_looking:
//		Compute L col by col
//		This is faster when it is memory bound.


namespace PackedCSparse {
	template <typename Tx, typename Ti>
	struct CholOutput : SparseMatrix<Tx, Ti>
	{
		TransposeOutput<bool, Ti> Lt;	// sparsity pattern of the Lt
		UniquePtr<Ti> diag;				// the index for diagonal element. Ax[diag[k]] is A_kk
		UniquePtr<Ti> c;				// c[i] = index the last nonzero on column i in the current L
		UniqueAlignedPtr<Tx> w;			// the row of L we are computing

		// The cost of this is roughly 3 times larger than chol
		// One can optimize it by using other data structure
		void initialize(const SparseMatrix<Tx, Ti>& A)
		{
			pcs_assert(A.initialized(), "chol: bad inputs.");
			pcs_assert(A.n == A.m, "chol: dimensions mismatch.");

			Ti n = A.n, * Ap = A.p.get(), * Ai = A.i.get();

			// initialize
			this->diag.reset(new Ti[n]);
			this->c.reset(new Ti[n]);
			this->w.reset(pcs_aligned_new<Tx>(n));

			// compute the sparsity pattern of L and diag
			using queue = std::priority_queue<Ti, std::vector<Ti>, std::greater<Ti>>;
			queue q; // sol'n of the current row of L
			Ti* mark = new Ti[n]; // used to prevent same indices push to q twice
			std::vector<Ti>* cols = new std::vector<Ti>[n]; // stores the non-zeros of each col of L
			Ti nz = 0, Anz = Ap[n];

			for (Ti i = 0; i < n; i++)
				mark[i] = -1;

			// for each row of A
			for (Ti i = 0; i < n; i++)
			{	// push i-th row of A, called a_12, into mark
				Ti s;
				for (s = Ap[i]; s < Ap[i + 1]; s++)
				{
					Ti j = Ai[s];
					if (j >= i) break;

					q.push(j);
					mark[j] = i;
				}
				if (s >= Anz) // this case happens only if the diagonal is 0. No cholesky in this case.
					this->diag[i] = 0;
				else
					this->diag[i] = s;

				// Solve L_11 l_12 = a_12
				while (!q.empty())
				{
					Ti j = q.top();

					for (Ti k : cols[j])
					{
						if (mark[k] != i)
						{
							q.push(k);
							mark[k] = i;
						}
					}
					q.pop();

					// update j col
					cols[j].push_back(i);
					++nz;
				}

				// diag
				cols[i].push_back(i);
				++nz;
			}
			delete[] mark;

			// write it as the compress form
			SparseMatrix<Tx, Ti>::initialize(n, n, nz);

			Ti s_start = 0; Ti s = 0;
			for (Ti i = 0; i < n; i++)
			{
				this->p[i] = s_start;
				for (Ti k : cols[i])
					this->i[s++] = k;
				s_start += Ti(cols[i].size());
			}
			this->p[n] = s_start;
			delete[] cols;

			this->Lt = transpose<Tx, Ti, bool>(*this);

			// initialize w to 0
			Tx Tv0 = Tx(0);
			for (Ti k = 0; k < n; k++)
				w[k] = Tv0;
		}
	};

	template <typename Tx, typename Ti>
	void chol(CholOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& A)
	{
		if (!o.initialized())
			o.initialize(A);

		//chol_up_looking(o, A);
		chol_left_looking(o, A);
	}

	template <typename Tx, typename Ti>
	void chol_up_looking(CholOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& A)
	{
		Ti *Ap = A.p.get(), * Ai = A.i.get(); Tx* Ax = A.x.get();
		Ti nzmax = o.nzmax; Ti n = A.n;
		Ti *Lp = o.p.get(); Ti* Li = o.i.get();
		Ti *Ltp = o.Lt.p.get(); Ti* Lti = o.Lt.i.get();

		Tx T0 = Tx(0);
		Tx* Lx = o.x.get(); Tx* w = o.w.get(); Ti* c = o.c.get();
		Ti* diag = o.diag.get();

		Ti* Lti_ptr = Lti;
		for (Ti k = 0; k < n; ++k)
		{
			c[k] = Lp[k];

			Ti s_end = diag[k];
			for (Ti s = Ap[k]; s < s_end; ++s)
				w[Ai[s]] = Ax[s];

			// Solve L_11 l_12 = a_12
			Tx d = Ax[s_end]; Ti i;
			for (; (i = *(Lti_ptr++)) < k;)
			{
				Ti dLi = Lp[i], ci = c[i]++;
				Tx Lki = w[i] / Lx[dLi];
				w[i] = T0; // maintain x = 0 for the (k+1) iteration

				for (Ti q = dLi + 1; q < ci; ++q)
					fnmadd(w[Li[q]], Lx[q], Lki);

				d -= Lki * Lki;
				Lx[ci] = Lki;
			}

			// l_22 = sqrt(a22 - <l12,l12>)
			Lx[c[k]++] = clipped_sqrt(d);
		}
	}

	template <typename Tx, typename Ti>
	void chol_left_looking(CholOutput<Tx, Ti>& o, const SparseMatrix<Tx, Ti>& A)
	{
		Ti* Ap = A.p.get(), * Ai = A.i.get(); Tx* Ax = A.x.get();
		Ti nzmax = o.nnz(); Ti n = A.n;
		Ti* Lp = o.p.get(); Ti* Li = o.i.get();
		Ti* Ltp = o.Lt.p.get(); Ti* Lti = o.Lt.i.get();

		Tx T0 = Tx(0), T1 = Tx(1);
		Tx* Lx = o.x.get();
		Tx* w = o.w.get(); Ti* c = o.c.get();
		Ti* diag = o.diag.get();

		for (Ti j = 0; j < n; ++j)
		{
			c[j] = Lp[j];

			// x = A_{j:n, j}
			{
				Ti is_start = diag[j], is_end = Ap[j + 1];
				for (Ti is = is_start; is < is_end; ++is)
					w[Ai[is]] = Ax[is];
			}

			// for each p in L_{j, 1:j-1}
			Ti ps_start = Ltp[j], ps_end = Ltp[j + 1] - 1;
			for (Ti ps = ps_start; ps < ps_end; ++ps)
			{
				Ti p = Lti[ps];
				Ti cp = c[p]++;
				Tx Ljp = Lx[cp];

				// for each i in L_{j:n,p}
				Ti is_start = cp, is_end = Lp[p + 1];
				for (Ti is = is_start; is < is_end; ++is)
				{
					Ti i = Li[is];
					fnmadd(w[i], Lx[is], Ljp);
				}
			}

			Tx Ljj = clipped_sqrt(w[j], 1e128);
			Lx[c[j]++] = Ljj;
			Tx inv_Ljj = T1 / Ljj;
			w[j] = T0;

			// for each i in L_{:,j}
			{
				Ti is_start = Lp[j] + 1, is_end = Lp[j + 1];
				for (Ti is = is_start; is < is_end; ++is)
				{
					Ti i = Li[is];
					Lx[is] = w[i] * inv_Ljj;
					w[i] = T0;
				}
			}
		}
	}

	template <typename Tx, typename Ti>
	CholOutput<Tx, Ti> chol(const SparseMatrix<Tx, Ti>& A)
	{
		CholOutput<Tx, Ti> o;
		chol(o, A);
		return o;
	}
}
