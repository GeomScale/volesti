/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLSMat.h.
   Arpack++ class ARluSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arlspen.h"

#ifndef ARLSMAT_H
#define ARLSMAT_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "blas1c.h"
#include "superluc.h"
#include "arlspdef.h"
#include "arlutil.h"

template<class ARTYPE> class ARluSymPencil;

template<class ARTYPE>
class ARluSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARluSymPencil<ARTYPE>;

 protected:

  bool        factored;
  char        uplo;
  int         order;
  int         nnz;
  int*        irow;
  int*        pcol;
  int*        permc;
  int*        permr;
  double      threshold;
  ARTYPE*     a;
  SuperMatrix A;
  SuperMatrix L;
  SuperMatrix U;
  ARhbMatrix<int, ARTYPE> mat;
  SuperLUStat_t stat;

  bool DataOK();

  virtual void Copy(const ARluSymMatrix& other);

  void ClearMem();

  void ExpandA(NCformat& A, NCformat& Aexp, ARTYPE sigma = (ARTYPE)0);

 public:

  int nzeros() { return nnz; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                    char uplop = 'L', double thresholdp = 0.1,
                    int orderp = 2, bool check = true);

  ARluSymMatrix();
  // Short constructor that does nothing.

  ARluSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
                char uplop = 'L', double thresholdp = 0.1,
                int orderp = 2, bool check = true);
  // Long constructor.

  ARluSymMatrix(const std::string& name, double thresholdp = 0.1,
                int orderp = 2, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARluSymMatrix(const ARluSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymMatrix() { ClearMem(); }
  // Destructor.

  ARluSymMatrix& operator=(const ARluSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARluSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
bool ARluSymMatrix<ARTYPE>::DataOK()
{

  int i, j, k;

  // Checking if pcol is in ascending order.

  i = 0;
  while ((i!=this->n)&&(pcol[i]<=pcol[i+1])) i++;
  if (i!=this->n) return false;

  // Checking if irow components are in order and within bounds.

  for (i=0; i!=this->n; i++) {
    j = pcol[i];
    k = pcol[i+1]-1;
    if (j<=k) {
      if (uplo == 'U') {
        if ((irow[j]<0)||(irow[k]>i)) return false;
      }
      else { // uplo == 'L'.
        if ((irow[j]<i)||(irow[k]>=this->n)) return false;
      }
      while ((j!=k)&&(irow[j]<irow[j+1])) j++;
      if (j!=k) return false;
    }
  }   

  return true;

} // DataOK.


template<class ARTYPE>
inline void ARluSymMatrix<ARTYPE>::Copy(const ARluSymMatrix<ARTYPE>& other)
{

  // Copying very fundamental variables.

  this->defined   = other.defined;
  factored  = other.factored;

  // Returning from here if "other" was not initialized.

  if (!this->defined) return;

  // Copying user-defined parameters.

  DefineMatrix(other.n, other.nnz, other.a, other.irow, other.pcol,
               other.uplo, other.threshold, other.order);

  // Throwing the original factorization away (this procedure 
  // is really awkward, but it is necessary because there
  // is no copy function for matrices L and U in the SuperLU 
  // library and it is not a good idea to do this kind of deep 
  // copy here).

  if (factored) {
    ArpackError(ArpackError::LAPACK_ERROR, "ARluSymMatrix");
    factored = false;
  }

} // Copy.


template<class ARTYPE>
void ARluSymMatrix<ARTYPE>::ClearMem()
{

  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
  }
  if (this->defined) {
    Destroy_SuperMatrix_Store(&A); // delete A.Store;
    delete[] permc;
    delete[] permr;
    permc = NULL;
    permr = NULL;
    A.Store = NULL;
  }

} // ClearMem.


template<class ARTYPE>
void ARluSymMatrix<ARTYPE>::
ExpandA(NCformat& A, NCformat& Aexp, ARTYPE sigma)
{

  // Defining local variables.

  bool   subtract;
  int    i, j, k;
  int    *colA, *colE;
  int    *indA, *indE;
  ARTYPE *valA, *valE;

  // Checking if sigma is zero.

  subtract = (sigma != (ARTYPE)0);

  // Simplifying the notation.

  valA = (ARTYPE*)A.nzval;
  valE = (ARTYPE*)Aexp.nzval;
  indA = (int*)A.rowind;
  indE = (int*)Aexp.rowind;
  colA = (int*)A.colptr;
  colE = (int*)Aexp.colptr;

  // Filling colE with zeros.

  for (i=0; i<=this->n; i++) colE[i] = 0;

  // Counting the elements in each column of A.

  if (uplo == 'U') {

    for (i=0; i!=this->n; i++) {
      k = colA[i+1];
      if ((k!=colA[i])&&(indA[k-1]==i)) {
        k--;
      }
      else {
        if (subtract) colE[i]++;
      }
      for (j=colA[i]; j<k; j++) colE[indA[j]]++;        
    }

  }
  else { // uplo == 'L'

    for (i=0; i!=this->n; i++) {
      k = colA[i];
      if ((k!=colA[i+1])&&(indA[k]==i)) {
        k++;
      }
      else {
        if (subtract) colE[i]++;
      }
      for (j=k; j<colA[i+1]; j++) colE[indA[j]]++;        
    }

  }  

  // Summing up colE elements.

  for (i=0; i<this->n; i++) colE[i+1]+=colE[i];

  // Adding colA to colE.

  for (i=this->n; i>0; i--) colE[i] = colE[i-1]+colA[i];
  colE[0] = colA[0];    

  // Expanding A.

  if (uplo == 'U') {

    for (i=0; i<this->n; i++) {
      for (j=colA[i]; j<(colA[i+1]-1); j++) {
        indE[colE[i]] = indA[j];
        indE[colE[indA[j]]] = i; 
        valE[colE[i]++] = valA[j];
        valE[colE[indA[j]]++] = valA[j];
      }
      if ((colA[i]!=colA[i+1])&&(indA[j]==i)) {
        indE[colE[i]] = i;
        if (subtract) {
          valE[colE[i]++] = valA[j]-sigma;
        }
        else {
          valE[colE[i]++] = valA[j];
        }
      }
      else {
        if (subtract) {
          indE[colE[i]] = i;
          valE[colE[i]++] = -sigma;
        }
      }
    }

  }
  else { // uplo  == 'L'

    for (i=0; i<this->n; i++) {
      k=colA[i];
      if ((k!=colA[i+1])&&(indA[k]==i)) {
        indE[colE[i]] = i;
        if (subtract) {
          valE[colE[i]++] = valA[k]-sigma;
        }
        else {
          valE[colE[i]++] = valA[k];
        }
        k++;
      }
      else {
        if (subtract) {
          indE[colE[i]] = i;
          valE[colE[i]++]  = -sigma;
        }
      }
      for (j=k; j<colA[i+1]; j++) {
        indE[colE[i]] = indA[j];
        indE[colE[indA[j]]] = i; 
        valE[colE[i]++] = valA[j];
        valE[colE[indA[j]]++] = valA[j];
      }
    }

  }

  // Adjusting index.

  for (i=this->n; i>0; i--) {
    colE[i] = colE[i-1];
  } 
  colE[0] = 0;

  Aexp.nnz = colE[this->n];

} // ExpandA.


template<class ARTYPE>
void ARluSymMatrix<ARTYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluSymMatrix::FactorA");
  }

  // Defining local variables.

  int         info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  ARTYPE*     aexp;
  SuperMatrix Aexp;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   Aexpstore;

  // Deleting previous versions of L and U.
  
  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
  }

  // Setting default values for gstrf parameters.

  int    panel_size = sp_ienv(1);
  int    relax      = sp_ienv(2);
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
  // options.DiagPivotThresh = 0.001;
  options.DiagPivotThresh = threshold;

  // Creating a temporary matrix Aexp.

  irowi = (int*)SUPERLU_MALLOC(sizeof(int) * (nnz*2));
  pcoli = (int*)SUPERLU_MALLOC(sizeof(int) * (this->n+1));
  aexp  = (ARTYPE*)SUPERLU_MALLOC(sizeof(ARTYPE) * (nnz*2));
  Create_CompCol_Matrix(&Aexp, this->n,  this->n, nnz, aexp, irowi, pcoli, SLU_NC, SLU_GE);

  // Expanding A.

  Astore    = (NCformat*)A.Store;
  Aexpstore = (NCformat*)Aexp.Store;
  ExpandA(*Astore, *Aexpstore);

  // Reserving memory for etree (used in matrix decomposition).

  etree = new int[this->n];

  // Defining LUStat.

  //StatInit(panel_size, relax);
  StatInit(&stat);

  // Defining the column permutation of matrix A
  // (using minimum degree ordering).

  get_perm_c(order, &Aexp, permc);

  // Permuting columns of A and creating the elimination tree.

  //sp_preorder("N", &Aexp, permc, etree, &AC);
  sp_preorder(&options, &Aexp, permc, etree, &AC);

  // Decomposing A.

//  gstrf("N",&AC, threshold, drop_tol, relax, panel_size, etree,
//        NULL, 0, permr, permc, &L, &U, &info);
  gstrf(&options,&AC, relax, panel_size, etree,
        NULL, 0, permc, permr, &L, &U, &stat, &info);

  // Deleting AC, Aexp and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&Aexp);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {        // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluSymMatrix::FactorA");
  }
  else if (info > this->n) {    // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluSymMatrix::FactorA");
  }
  else if (info > 0) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluSymMatrix::FactorA");
  }

} // FactorA.


template<class ARTYPE>
void ARluSymMatrix<ARTYPE>::FactorAsI(ARTYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluSymMatrix::FactorAsI");
  }

  // Defining local variables.

  int         info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  ARTYPE*     asi;
  SuperMatrix AsI;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   AsIstore;

  // Deleting previous versions of L and U.
  
  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
  }

  // Setting default values for gstrf parameters.

  int    panel_size = sp_ienv(1);
  int    relax      = sp_ienv(2);
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
  // options.DiagPivotThresh = 0.001;
  options.DiagPivotThresh = threshold;

  // Creating a temporary matrix AsI.

  irowi = (int*)SUPERLU_MALLOC(sizeof(int) * (nnz*2+this->n));
  pcoli = (int*)SUPERLU_MALLOC(sizeof(int) * (this->n+1));
  asi   = (ARTYPE*)SUPERLU_MALLOC(sizeof(ARTYPE) * (nnz*2+this->n));
  Create_CompCol_Matrix(&AsI, this->n,  this->n, nnz, asi, irowi, pcoli, SLU_NC, SLU_GE);

  // Subtracting sigma*I from A and storing the result on AsI.

  Astore   = (NCformat*)A.Store;
  AsIstore = (NCformat*)AsI.Store;
  ExpandA(*Astore, *AsIstore, sigma);

  // Reserving memory for etree (used in matrix decomposition).

  etree = new int[this->n];

  // Defining LUStat.

  //StatInit(panel_size, relax);
  StatInit(&stat);

  // Defining the column permutation of matrix AsI
  // (using minimum degree ordering).

  get_perm_c(order, &AsI, permc);

  // Permuting columns of AsI and creating the elimination tree.

  sp_preorder(&options, &AsI, permc, etree, &AC);

  // Decomposing AsI.

//  gstrf("N",&AC, threshold, drop_tol, relax, panel_size, etree,
//        NULL, 0, permr, permc, &L, &U, &info);
  gstrf(&options,&AC, relax, panel_size, etree,
        NULL, 0, permc, permr, &L, &U, &stat, &info);

  // Deleting AC, AsI and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&AsI);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {        // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluSymMatrix::FactorAsI");
  }
  else if (info > this->n) {    // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluSymMatrix::FactorAsI");
  }
  else if (info > 0) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluSymMatrix::FactorAsI");
  }

} // FactorAsI.


template<class ARTYPE>
void ARluSymMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{

  int    i, j, k;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluSymMatrix::MultMv");
  }

  // Determining w = M.v.

  for (i=0; i!=this->m; i++) w[i]=(ARTYPE)0;

  if (uplo == 'U') {

    for (i=0; i!=this->n; i++) {
      t = v[i];
      k = pcol[i+1];
      if ((k!=pcol[i])&&(irow[k-1]==i)) {
        w[i] += t*a[k-1];
        k--;
      }
      for (j=pcol[i]; j<k; j++) {
        w[irow[j]] += t*a[j];
        w[i] += v[irow[j]]*a[j];
      }
    }

  }
  else {

    for (i=0; i!=this->n; i++) {
      t = v[i];
      k = pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        w[i] += t*a[k];
        k++;
      }
      for (j=k; j<pcol[i+1]; j++) {
        w[irow[j]] += t*a[j];
        w[i] += v[irow[j]]*a[j];
      }
    }

  }

} // MultMv.


template<class ARTYPE>
void ARluSymMatrix<ARTYPE>::MultInvv(ARTYPE* v, ARTYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARluSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).

  int         info;
  SuperMatrix B;

  if (&v != &w) copy(this->n, v, 1, w, 1);
  Create_Dense_Matrix(&B, this->n, 1, w, this->n, SLU_DN, SLU_GE);
//  gstrs("N", &L, &U, permr, permc, &B, &info);
  StatInit(&stat);
  trans_t trans = NOTRANS;
  gstrs(trans, &L, &U, permc, permr, &B, &stat, &info);
  Destroy_SuperMatrix_Store(&B); // delete B.Store;

} // MultInvv.


template<class ARTYPE>
inline void ARluSymMatrix<ARTYPE>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp, int* pcolp,
             char uplop, double thresholdp, int orderp, bool check)
{

  this->m         = np;
  this->n         = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[this->n]   = nnz;
  uplo      = uplop;
  threshold = thresholdp;
  order     = orderp;

  // Checking data.

  if ((check)&&(!DataOK())) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARluSymMatrix::DefineMatrix");
  }

  // Creating SuperMatrix A.

  Create_CompCol_Matrix(&A, this->n, this->n, nnz, a, irow, pcol, SLU_NC, SLU_GE);

  // Reserving memory for vectors used in matrix decomposition.

  permc = new int[this->n];
  permr = new int[this->n];

  this->defined = true;

} // DefineMatrix.


template<class ARTYPE>
inline ARluSymMatrix<ARTYPE>::ARluSymMatrix(): ARMatrix<ARTYPE>()
{

  factored = false;
  permc    = NULL;
  permr    = NULL;
 
} // Short constructor.


template<class ARTYPE>
inline ARluSymMatrix<ARTYPE>::
ARluSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
              int* pcolp, char uplop, double thresholdp,
              int orderp, bool check)                   : ARMatrix<ARTYPE>(np)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, uplop, thresholdp, orderp, check);

} // Long constructor.


template<class ARTYPE>
ARluSymMatrix<ARTYPE>::
ARluSymMatrix(const std::string& file, double thresholdp, int orderp, bool check)
{

  factored = false;

  try {
    mat.Define(file);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARluSymMatrix");
  }

  if ((mat.NCols() == mat.NRows()) && (mat.IsSymmetric())) {

    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), 'L', thresholdp, orderp, check);
  }
  else {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARluSymMatrix::ARluSymMatrix");
  }

} // Long constructor (Harwell-Boeing file).


template<class ARTYPE>
ARluSymMatrix<ARTYPE>& ARluSymMatrix<ARTYPE>::
operator=(const ARluSymMatrix<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSMAT_H
