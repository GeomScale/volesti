/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE UMFPACKc.h.
   Interface to UMFPACK routines.

   Author of this class:
      Martin Reuter
      Date 2/28/2013
      
   Arpack++ Author:
      Francisco Gomes
      
   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef UMFPACKC_H
#define UMFPACKC_H

#ifdef __cplusplus
extern "C" {
#endif

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20
#define UMFPACK_OK (0)
#define UMFPACK_A	(0)	/* Ax=b    */
#define UMFPACK_PRL 0			/* print level */

void umfpack_di_defaults
(
    double Control [UMFPACK_CONTROL]
) ;


int umfpack_di_symbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_di_numeric
(
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

void umfpack_di_free_symbolic
(
    void **Symbolic
) ;

void umfpack_di_free_numeric
(
    void **Numeric
) ;

int umfpack_di_triplet_to_col
(
    int n_row,
    int n_col,
    int nz,
    const int Ti [ ],
    const int Tj [ ],
    const double Tx [ ],
    int Ap [ ],
    int Ai [ ],
    double Ax [ ],
    int Map [ ]
) ;

int umfpack_di_solve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_di_report_matrix
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    int col_form,
    const double Control [UMFPACK_CONTROL]
) ;

#ifdef __cplusplus
  }
#endif

//#include "umfpack.h"
#include <fstream>

inline void Write_Triplet_Matrix(const std::string & fname, int * tripi,
                                 int * tripj, double* tripx, unsigned int nnz)
{
  std::ofstream myfile; 
  myfile.open ( fname.c_str() );
	myfile.precision(20);
  for (unsigned int i=0;i<nnz;i++)
  {
    myfile << tripi[i]+1 << " " << tripj[i]+1 << " " << tripx[i] << std::endl;
  }
  myfile.close();
}

/*inline void Write_Cholmod_Sparse_Matrix(const std::string & fname,
                             cholmod_sparse* A, cholmod_common *c)
{
  std::ofstream myfile; 
  myfile.open ( fname.c_str() );
  cholmod_triplet * T = cholmod_sparse_to_triplet(A,c);
  //std::cout << " [ " << std::endl;
	myfile.precision(20);
  for (unsigned int i=0;i<T->nnz;i++)
  {
    myfile << ((int*)T->i)[i]+1 << " " << ((int*)T->j)[i]+1 << " " << ((double*)T->x)[i] << std::endl;
  }
  //std::cout << " ] " << std::endl;
  myfile.close();
  
  cholmod_free_triplet(&T,c);

}

// Create_Cholmod_Sparse_Matrix 
inline cholmod_sparse* Create_Cholmod_Sparse_Matrix(int m, int n, int nnz,
      double* a, int* irow, int* pcol, char uplo, cholmod_common *c)
{
  
  cholmod_sparse* A = new cholmod_sparse;
  A->nrow = m;
  A->ncol = n;
  A->nzmax = nnz;
  A->p = pcol;
  A->i = irow;
  A->nz = NULL;
  A->x = a;
  A->z = NULL;
  if (uplo == 'L') A->stype = -1;
  else A->stype = 1;
  A->itype = CHOLMOD_INT;
  A->xtype = CHOLMOD_REAL; // real
  A->dtype = CHOLMOD_DOUBLE; // double
  A->sorted = 0;
  A->packed = 1;

  return A;  
  
  

  
} // Create_Cholmod_Sparse_Matrix (double).

// Create_Cholmod_Dense_Matrix (from Triplet)
inline cholmod_dense* Create_Cholmod_Dense_Matrix(int m, int n,
                                  double* a, cholmod_common *c)
{


  cholmod_dense* A = new cholmod_dense;
  A->nrow = m;
  A->ncol = n;
  A->nzmax = m*n;
  A->d = m;
  A->x = a;
  A->z = NULL;
  A->xtype = CHOLMOD_REAL; // real
  A->dtype = CHOLMOD_DOUBLE; // double

//  cholmod_dense* As = cholmod_copy_dense(A,c);
  
  return A;
  
} // Create_Cholmod_Dense_Matrix (double).

// Create_Cholmod_Dense_Matrix (from Triplet)
inline void Get_Cholmod_Dense_Data(cholmod_dense* A, int n, double* a)
{
  memcpy(a,A->x,n*sizeof(double));
  
//  for (int i = 0;i<n;i++)
//    a[i] = ((double*)A->x)[i];
  
} // Create_Cholmod_Dense_Matrix (double).

*/

#endif // UMFPACKC_H
