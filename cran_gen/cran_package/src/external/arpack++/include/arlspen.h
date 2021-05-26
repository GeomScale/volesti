/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLSPen.h.
   Arpack++ class ARluSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLSPEN_H
#define ARLSPEN_H

#include <cstddef>

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "superluc.h"
#include "arlspdef.h"
#include "arlutil.h"
#include "arlsmat.h"


template<class ARTYPE>
class ARluSymPencil
{

 protected:

  bool                   factored;
  int*                   permc;
  int*                   permr;
  char                   part;
  char                   uplo;
  ARluSymMatrix<ARTYPE>* A;
  ARluSymMatrix<ARTYPE>* B;
  SuperMatrix            L;
  SuperMatrix            U;
  SuperLUStat_t stat;

  virtual void Copy(const ARluSymPencil& other);

  void ClearMem();

  void SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
                   int yind[], int ny, ARTYPE z[], int zind[], int& nz);

  void ExpandAsB(int n, NCformat& AsB);

  void SubtractAsB(int n, ARTYPE sigma, NCformat& A, 
                   NCformat& B, NCformat& AsB);

 public:

  bool IsFactored() { return factored; }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w);

  void DefineMatrices(ARluSymMatrix<ARTYPE>& Ap, ARluSymMatrix<ARTYPE>& Bp);

  ARluSymPencil();
  // Short constructor that does nothing.

  ARluSymPencil(ARluSymMatrix<ARTYPE>& Ap, ARluSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARluSymPencil(const ARluSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymPencil() { ClearMem(); }
  // Destructor.

  ARluSymPencil& operator=(const ARluSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARluSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARluSymPencil<ARTYPE>::
Copy(const ARluSymPencil<ARTYPE>& other)
{

  factored = other.factored;
  part     = other.part;
  uplo     = other.uplo;
  A        = other.A;
  B        = other.B;

  // Throwing the original factorization away (this procedure 
  // is really awkward, but it is necessary because there
  // is no copy function for matrices L and U in the SuperLU 
  // library and it is not a good idea to do this kind of deep 
  // copy here).

  if (factored) {
    ArpackError(ArpackError::DISCARDING_FACTORS, "ARluSymPencil");
    factored = false;
  }

} // Copy.


template<class ARTYPE>
void ARluSymPencil<ARTYPE>::ClearMem()
{

  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
    delete[] permc;
    delete[] permr;
    permc = NULL;
    permr = NULL;
  }

} // ClearMem.


template<class ARTYPE>
void ARluSymPencil<ARTYPE>::
SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
            int yind[], int ny, ARTYPE z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  int ix, iy;

  nz = 0;
  if ((nx == 0) || (a == (ARTYPE)0)) {
    copy(ny,y,1,z,1);
    for (iy=0; iy!=ny; iy++) zind[iy] = yind[iy];
    nz = ny;
    return;
  }
  if (ny == 0) {
    copy(nx,x,1,z,1);
    scal(nx,a,z,1);
    for (ix=0; ix!=nx; ix++) zind[ix] = xind[ix];
    nz = nx;
    return;
  }
  ix = 0;
  iy = 0;
  while (true) {
    if (xind[ix] == yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++]+y[iy++];
      if ((ix == nx)||(iy == ny)) break;
    }
    else if (xind[ix] < yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++];
      if (ix == nx) break;
    }
    else {
      zind[nz] = yind[iy];
      z[nz++]  = y[iy++];
      if (iy == ny) break;
    }
  }
  while (iy < ny) {
    zind[nz] = yind[iy];
    z[nz++]  = y[iy++];
  }
  while (ix < nx) {
    zind[nz] = xind[ix];
    z[nz++]  = x[ix++];
  }

} // SparseSaxpy.


template<class ARTYPE>
void ARluSymPencil<ARTYPE>::ExpandAsB(int n, NCformat& AsB)
{

  int    i, j, k;
  int    *pcol, *pos, *col, *ind;
  ARTYPE *val;

  // simplifying the notation.
  
  val = (ARTYPE*)AsB.nzval;
  ind = AsB.rowind;
  col = AsB.colptr;

  // Initializing vectors.

  pcol   = new int[n+1];
  pos    = new int[n+1];
  for (i=0; i<=n; i++) pcol[i] = col[i];
  for (i=0; i<=n; i++) pos[i]  = 0;

  // Counting the elements in each column of AsB.

  if (uplo == 'U') {

    for (i=0; i!=n; i++) {
      k = pcol[i+1];
      if ((k!=pcol[i])&&(ind[k-1]==i)) k--;
      for (j=pcol[i]; j<k; j++) pos[ind[j]]++;        
    }

  }
  else { // uplo == 'L'

    for (i=0; i!=n; i++) {
      k = pcol[i];
      if ((k!=pcol[i+1])&&(ind[k]==i)) k++;
      for (j=k; j<pcol[i+1]; j++) pos[ind[j]]++;        
    }

  }  

  // Summing up index elements.

  for (i=0; i<n; i++) pos[i+1] += pos[i];
  for (i=n; i>0; i--) col[i] += pos[i-1];
    
  // Expanding A.

  if (uplo == 'U') {

    for (i=n-1; i>=0; i--) {
      pos[i] = col[i]+pcol[i+1]-pcol[i];
      k = pos[i]-1;
      for (j=pcol[i+1]-1; j>=pcol[i]; j--) {
        val[k]   = val[j];
        ind[k--] = ind[j];
      }
    }
    for (i=1; i<n; i++) {
      k = col[i]+pcol[i+1]-pcol[i];
      if ((k>col[i])&&(ind[k-1]==i)) k--;
      for (j=col[i]; j<k; j++) {
        val[pos[ind[j]]]   = val[j];
        ind[pos[ind[j]]++] = i;
      }
    }

  }
  else { // uplo  == 'L'

    for (i=n-1; i>=0; i--) {
      k = col[i+1]-1;
      for (j=pcol[i+1]-1; j>=pcol[i]; j--) {
        val[k]   = val[j];
        ind[k--] = ind[j];
      }
      pos[i] = col[i];
    }
    for (i=0; i<(n-1); i++) {
      k = col[i+1]-pcol[i+1]+pcol[i];
      if ((k<col[i+1])&&(ind[k]==i)) k++;
      for (j=k; j<col[i+1]; j++) {
        val[pos[ind[j]]]   = val[j];
        ind[pos[ind[j]]++] = i;
      }
    }

  }

  AsB.nnz = col[n]; 

  //  Deleting temporary vectors.

  delete[] pcol;  
  delete[] pos;

} // ExpandAsB.


template<class ARTYPE>
void ARluSymPencil<ARTYPE>::
SubtractAsB(int n, ARTYPE sigma, NCformat& matA, NCformat& matB, NCformat& AsB)
{

  int     i, acol, bcol, asbcol, scol;
  ARTYPE* anzval;
  ARTYPE* bnzval;
  ARTYPE* asbnzval;

  // Quitting function if A->uplo is not equal to B->uplo.

  if ((A->uplo != B->uplo)&&(sigma != (ARTYPE)0)) {
    throw ArpackError(ArpackError::DIFFERENT_TRIANGLES,
                      "ARluSymPencil::SubtractAsB");
  }
  uplo = A->uplo;

  // Telling the compiler that nzval must ve viewed as a vector of ARTYPE.

  anzval   = (ARTYPE*)matA.nzval;
  bnzval   = (ARTYPE*)matB.nzval;
  asbnzval = (ARTYPE*)AsB.nzval;

  // Subtracting sigma*B from A.

  AsB.colptr[0] = 0;
  asbcol = 0;

  for (i=0; i!=n; i++) {
    bcol = matB.colptr[i];
    acol = matA.colptr[i];
    SparseSaxpy(-sigma, &bnzval[bcol], &matB.rowind[bcol], 
                matB.colptr[i+1]-bcol, &anzval[acol], &matA.rowind[acol], 
                matA.colptr[i+1]-acol, &asbnzval[asbcol], 
                &AsB.rowind[asbcol], scol);
    asbcol += scol;
    AsB.colptr[i+1] = asbcol;
  }

  // Expanding AsB.

  ExpandAsB(n, AsB);

} // SubtractAsB.


// template<class ARTYPE>
// void ARluSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
// {
// 
//   // Quitting the function if A and B were not defined.
// 
//   if (!(A->IsDefined()&&B->IsDefined())) {
//     throw ArpackError(ArpackError::DATA_UNDEFINED,
//                       "ARluSymPencil::FactorAsB");
//   }
// 
//   // Quitting the function if A and B are not square.
// 
//   if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
//     throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
//                       "ARluSymPencil::FactorAsB");
//   }
// 
//   // Defining local variables.
// 
//   int         nnzi, info;
//   int*        etree;
//   int*        irowi;
//   int*        pcoli;
//   ARTYPE*     asb;
//   SuperMatrix AsB;
//   SuperMatrix AC;
//   NCformat*   Astore;
//   NCformat*   Bstore;
//   NCformat*   AsBstore;
// 
//   // Deleting old versions of L, U, perm_r and perm_c.
// 
//   ClearMem();
// 
//   // Setting default values for gstrf parameters.
// 
//   ARTYPE drop_tol = (ARTYPE)0;
//   int  panel_size = sp_ienv(1);
//   int  relax      = sp_ienv(2);
// 
//   // Defining A and B format.
// 
//   Astore = (NCformat*)A->A.Store;
//   Bstore = (NCformat*)B->A.Store;
// 
//   // Creating a temporary matrix AsB.
// 
//   nnzi  = (Astore->nnz+Bstore->nnz)*2;
//   irowi = new int[nnzi];
//   pcoli = new int[A->ncols()+1];
//   asb   = new ARTYPE[nnzi];
//   Create_CompCol_Matrix(&AsB, A->nrows(), A->ncols(), nnzi, asb,
//                         irowi, pcoli, NC, GE);
// 
//   // Subtracting sigma*B from A and storing the result on AsB.
// 
//   AsBstore = (NCformat*)AsB.Store;
//   SubtractAsB(A->ncols(), sigma, *Astore, *Bstore, *AsBstore);
// 
//   // Reserving memory for some vectors used in matrix decomposition.
// 
//   etree = new int[A->ncols()];
//   if (permc == NULL) permc = new int[A->ncols()];
//   if (permr == NULL) permr = new int[A->ncols()];
// 
//   // Defining LUStat.
// 
//   StatInit(panel_size, relax);
// 
//   // Defining the column permutation of matrix AsB
//   // (using minimum degree ordering on AsB'*AsB).
// 
//   get_perm_c(A->order, &AsB, permc);
// 
//   // Permuting columns of AsB and
//   // creating the elimination tree of AsB'*AsB.
// 
//   sp_preorder("N", &AsB, permc, etree, &AC);
// 
//   // Decomposing AsB.
// 
//   gstrf("N",&AC, A->threshold, drop_tol, relax, panel_size, etree,
//         NULL, 0, permr, permc, &L, &U, &info);
// 
//   // Deleting AC, AsB and etree.
// 
//   Destroy_CompCol_Permuted(&AC);
//   Destroy_CompCol_Matrix(&AsB);
//   delete[] etree;
// 
//   factored = (info == 0);
// 
//   // Handling errors.
// 
//   if (info < 0)  {              // Illegal argument.
//     throw ArpackError(ArpackError::PARAMETER_ERROR,
//                       "ARluSymPencil::FactorAsB");
//   }
//   else if (info > A->ncols()) {  // Memory is not sufficient.
//     throw ArpackError(ArpackError::MEMORY_OVERFLOW,
//                       "ARluSymPencil::FactorAsB");
//   }
//   else if (info > 0) {          // Matrix is singular.
//     throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
//                       "ARluSymPencil::FactorAsB");
//   }
// 
// } // FactorAsB.

template<class ARTYPE>
void ARluSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARluSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARluSymPencil::FactorAsB");
  }

  // Defining local variables.

  int         nnzi, info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  ARTYPE*     asb;
  SuperMatrix AsB;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   Bstore;
  NCformat*   AsBstore;

  // Deleting old versions of L, U, perm_r and perm_c.

  ClearMem();

  // Setting default values for gstrf parameters.

  int  panel_size = sp_ienv(1);
  int  relax      = sp_ienv(2);
  superlu_options_t options;
  /* Set the default input options:
  options.Fact = DOFACT;
  options.Equil = YES;
  options.ColPerm = COLAMD;
  options.DiagPivotThresh = 1.0;
  options.Trans = NOTRANS;
  options.IterRefine = NOREFINE;
  options.SymmetricMode = NO;
  options.PivotGrowth = NO;
  options.ConditionNumber = NO;
  options.PrintStat = YES;
  */
  set_default_options(&options);

  /* Now we modify the default options to use the symmetric mode. */
  options.SymmetricMode = YES;
  options.ColPerm = MMD_AT_PLUS_A;
  options.DiagPivotThresh = A->threshold;

  // Defining A and B format.

  Astore = (NCformat*)A->A.Store;
  Bstore = (NCformat*)B->A.Store;

  // Creating a temporary matrix AsB.

  nnzi  = (Astore->nnz+Bstore->nnz)*2;
  irowi = new int[nnzi];
  pcoli = new int[A->ncols()+1];
  asb   = new ARTYPE[nnzi];
  Create_CompCol_Matrix(&AsB, A->nrows(), A->ncols(), nnzi, asb,
                        irowi, pcoli, SLU_NC, SLU_GE);

  // Subtracting sigma*B from A and storing the result on AsB.

  AsBstore = (NCformat*)AsB.Store;
  SubtractAsB(A->ncols(), sigma, *Astore, *Bstore, *AsBstore);

  // Reserving memory for some vectors used in matrix decomposition.

  etree = new int[A->ncols()];
  if (permc == NULL) permc = new int[A->ncols()];
  if (permr == NULL) permr = new int[A->ncols()];

  // Defining LUStat.

//  StatInit(panel_size, relax);
    SuperLUStat_t stat;
    StatInit(&stat);

  // Defining the column permutation of matrix AsB
  // (using minimum degree ordering on AsB'*AsB).

  get_perm_c(A->order, &AsB, permc);

  // Permuting columns of AsB and
  // creating the elimination tree of AsB'*AsB.

  sp_preorder(&options, &AsB, permc, etree, &AC);

  // Decomposing AsB.

//  gstrf("N",&AC, A->threshold, drop_tol, relax, panel_size, etree,
//        NULL, 0, permr, permc, &L, &U, &info);
  gstrf(&options, &AC, relax, panel_size, etree,
        NULL, 0, permc, permr, &L, &U, &stat, &info);

  // Deleting AC, AsB and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&AsB);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {              // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluSymPencil::FactorAsB");
  }
  else if (info > A->ncols()) {  // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluSymPencil::FactorAsB");
  }
  else if (info > 0) {          // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluSymPencil::FactorAsB");
  }

} // FactorAsB.

template<class ARTYPE>
void ARluSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{
  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.


template<class ARTYPE>
void ARluSymPencil<ARTYPE>::MultInvAsBv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if AsB was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARluSymPencil::MultInvAsBv");
  }

  // Solving AsB.w = v.

  int         info;
  SuperMatrix RHS;

  copy(A->nrows(), v, 1, w, 1);
  Create_Dense_Matrix(&RHS, A->nrows(), 1, w, A->nrows(), SLU_DN, SLU_GE);
//  gstrs("N", &L, &U, permr, permc, &RHS, &info);
  trans_t trans = NOTRANS;
  StatInit(&stat);

  gstrs(trans, &L, &U, permc, permr, &RHS, &stat, &info);

  Destroy_SuperMatrix_Store(&RHS); // delete RHS.Store;

} // MultInvAsBv.


template<class ARTYPE>
inline void ARluSymPencil<ARTYPE>::
DefineMatrices(ARluSymMatrix<ARTYPE>& Ap, ARluSymMatrix<ARTYPE>& Bp)
{

  A     = &Ap;
  B     = &Bp;
  permc = NULL;
  permr = NULL;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARluSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE>
inline ARluSymPencil<ARTYPE>::ARluSymPencil()
{
  
  factored = false; 
  part     = 'N'; 
  permr    = NULL;
  permc    = NULL;

} // Short constructor.


template<class ARTYPE>
inline ARluSymPencil<ARTYPE>::
ARluSymPencil(ARluSymMatrix<ARTYPE>& Ap, ARluSymMatrix<ARTYPE>& Bp)
{

  factored = false;
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE>
ARluSymPencil<ARTYPE>& ARluSymPencil<ARTYPE>::
operator=(const ARluSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSPEN_H
