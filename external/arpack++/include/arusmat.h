/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSMat.h.
   Arpack++ class ARumSymMatrix definition.

   Modified to work with Umfpack v5.??
      Martin Reuter
      Date 02/28/2013

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


#include "aruspen.h"

#ifndef ARUSMAT_H
#define ARUSMAT_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
//#include "blas1c.h"
#include "umfpackc.h"

template<class ARTYPE> class ARumSymPencil;

template<class ARTYPE>
class ARumSymMatrix: public ARMatrix<ARTYPE> {

  friend class ARumSymPencil<ARTYPE>;

 protected:

  bool    factored;
  char    uplo;
  int     nnz;
 /* int     fillin;
  int     lvalue;
  int     lindex;
  int     keep[20];
  int     icntl[20];
  int     info[40];
  ARTYPE  cntl[10];
  ARTYPE  rinfo[20];
  int*    index;
  ARTYPE* value;*/
  int*    irow;
  int*    pcol;
  int     status;
  double  threshold;
  ARTYPE* a;
  ARhbMatrix<int, ARTYPE> mat;
  void*   Numeric;
  int*    Ap;
  int*    Ai;
  ARTYPE* Ax; 

  bool DataOK();

  virtual void Copy(const ARumSymMatrix& other);

  void ClearMem();

  void ExpandA(ARTYPE sigma = (ARTYPE)0);

//  void CreateStructure();

  void ThrowError();

 public:

  int nzeros() { return nnz; }

//  int  FillFact() { return fillin; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(ARTYPE sigma);

  void MultMv(ARTYPE* v, ARTYPE* w);

  void MultInvv(ARTYPE* v, ARTYPE* w);

  void DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                    int* pcolp, char uplop = 'L', double thresholdp = 0.1, 
                    int fillinp = 9, bool reducible = true, bool check = true);

  ARumSymMatrix(): ARMatrix<ARTYPE>()
  {
    factored = false;
    Numeric = NULL;
    Ap = NULL;
    Ai = NULL;
    Ax = NULL;
  }
  // Short constructor that does nothing.

  ARumSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
                int* pcolp, char uplop = 'L', double thresholdp = 0.1,
                int fillinp = 9, bool reducible = true, bool check = true);
  // Long constructor.

  ARumSymMatrix(const std::string& name, double thresholdp = 0.1, int fillinp = 9,
                bool reducible = true, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARumSymMatrix(const ARumSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymMatrix() { ClearMem(); }
  // Destructor.

  ARumSymMatrix& operator=(const ARumSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
bool ARumSymMatrix<ARTYPE>::DataOK()
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
inline void ARumSymMatrix<ARTYPE>::ClearMem()
{

  if (factored)
  {
    if (Numeric) umfpack_di_free_numeric (&Numeric);
    //if (value) delete[] value;
    //if (index) delete[] index;
    //value = NULL;
    //index = NULL;
    if (Ai) delete [] Ai;
    Ai = NULL;
    if (Ap) delete [] Ap;
    Ap = NULL;
    if (Ax) delete [] Ax;
    Ax = NULL;
  }

} // ClearMem.



template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::Copy(const ARumSymMatrix<ARTYPE>& other)
{

  // Copying very fundamental variables.
  ClearMem();

  // Copying very fundamental variables and user-defined parameters.

  this->m         = other.m;
  this->n         = other.n;
  this->defined   = other.defined;
  factored  = other.factored;
  //fillin    = other.fillin;
  nnz       = other.nnz;
  //lvalue    = other.lvalue;
 //lindex    = other.lindex;
  irow      = other.irow;
  pcol      = other.pcol;
  a         = other.a;
  threshold = other.threshold;
  uplo      = other.uplo;

  // Returning from here if "other" was not initialized.

  if (!this->defined) return;

  // Returning from here if "other" was not factored.

  if (!factored) return;

  factored = false;

} // Copy.

template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::ExpandA(ARTYPE sigma)
{
std::cout <<"ARumSymMatrix::ExpandA(" << sigma << ") ..." << std::flush; 

  ClearMem();
 
  // Checking if sigma is zero.
  bool subtract = (sigma != (ARTYPE)0);

  int mynnz = 2*nnz;
  if (subtract) mynnz = 2*nnz + this->n; // some space for the diag entries just in case
  
  // create triples (i,j,value)
  int * tripi = new int[mynnz];
  int * tripj = new int[mynnz];
  ARTYPE* tripx = new ARTYPE[mynnz];
  int count = 0;
  int i,j;
//  if (uplo == 'U')
  {
    for (i=0; i != this->n; i++)
    {
      bool founddiag = false;
      for (j=pcol[i]; j<(pcol[i+1]); j++)
      {
        
        if (i == irow[j]) // on diag
        {
          tripi[count] = i;
          tripj[count] = irow[j];
          if (subtract)
          {
            tripx[count] = a[j]-sigma;
            founddiag = true;
          }
          else tripx[count] = a[j];
          count++;
        }
        else
        {
        
          tripi[count] = i;
          tripj[count] = irow[j];
          tripx[count] = a[j];
          count++;
          tripj[count] = i;
          tripi[count] = irow[j];
          tripx[count] = a[j];
          count++;
        }
      }
      if (subtract && ! founddiag)
      {
        tripi[count] = i;
        tripj[count] = i;
        tripx[count] = -sigma;
        count++;
      }
    }
  }
  
  // convert triples to Ax Ap Ai
  Ap = new int[this->n+1];
  Ai = new int[count];
  Ax = new ARTYPE[count];
  status = umfpack_di_triplet_to_col (this->n, this->n, count, tripi, tripj, tripx, Ap, Ai, Ax,  (int *)NULL) ;
  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymMatrix::ExpandA");
  if (Ap[this->n] != count)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymMatrix::ExpandA");


  // cleanup
  delete [] tripi;
  delete [] tripj;
  delete [] tripx;

  //std::cout << std::endl << std::endl;
  //double Control [UMFPACK_CONTROL];
  //Control [UMFPACK_PRL] = 3;
  //status = umfpack_di_report_matrix(this->n, this->n,Ap, Ai, Ax,0,Control);
  //std::cout << " status: " << status << std::endl;
  //std::cout << std::endl << std::endl;

  std::cout <<" done!" << std::endl; 

}

/*template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::ExpandA(ARTYPE sigma)
{

  bool subtract;
  int  i, j, k, ki;

  // Checking if sigma is zero.

  subtract = (sigma != (ARTYPE)0);

  // Filling index with zeros.

  for (i=0; i<=this->n; i++) index[i] = 0;

  // Counting the elements in each column of A.

  if (uplo == 'U') {

    for (i=0; i!=this->n; i++) {
      k = pcol[i+1];
      if ((k!=pcol[i])&&(irow[k-1]==i)) {
        k--;
      }
      else {
        if (subtract) index[i]++;
      }
      for (j=pcol[i]; j<k; j++) index[irow[j]]++;        
    }

  }
  else { // uplo == 'L'

    for (i=0; i!=this->n; i++) {
      k = pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        k++;
      }
      else {
        if (subtract) index[i]++;
      }
      for (j=k; j<pcol[i+1]; j++) index[irow[j]]++;        
    }

  }  

  // Summing up index elements.

  for (i=0; i<this->n; i++) index[i+1]+=index[i];

  // Adding pcol to index.

  for (i=this->n; i>0; i--) index[i] = index[i-1]+pcol[i];
  index[0] = pcol[0];    

  // Expanding A.

  ki = this->n+1;

  if (uplo == 'U') {

    for (i=0; i<this->n; i++) {
      for (j=pcol[i]; j<(pcol[i+1]-1); j++) {
        index[ki+index[i]] = irow[j]+1;
        index[ki+index[irow[j]]] = i+1; 
        value[index[i]++] = a[j];
        value[index[irow[j]]++] = a[j];
      }
      if ((pcol[i]!=pcol[i+1])&&(irow[j]==i)) {
        index[ki+index[i]] = i+1;
        if (subtract) {
          value[index[i]++] = a[j]-sigma;
        }
        else {
          value[index[i]++] = a[j];
        }
      }
      else {
        if (subtract) {
          index[ki+index[i]] = i+1;
          value[index[i]++]  = -sigma;
        }
      }
    }

  }
  else { // uplo  == 'L'

    for (i=0; i<this->n; i++) {
      k=pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        index[ki+index[i]] = i+1;
        if (subtract) {
          value[index[i]++] = a[k]-sigma;
        }
        else {
          value[index[i]++] = a[k];
        }
        k++;
      }
      else {
        if (subtract) {
          index[ki+index[i]] = i+1;
          value[index[i]++]  = -sigma;
        }
      }
      for (j=k; j<pcol[i+1]; j++) {
        index[ki+index[i]] = irow[j]+1;
        index[ki+index[irow[j]]] = i+1; 
        value[index[i]++] = a[j];
        value[index[irow[j]]++] = a[j];
      }
    }

  }

  // Adjusting index.

  for (i=this->n; i>0; i--) {
    index[i] = index[i-1]+1;
  } 
  index[0] = 1;

} // ExpandA.*/


/*template<class ARTYPE>
inline void ARumSymMatrix<ARTYPE>::CreateStructure()
{

  int dimfact = (((fillin+1)*nnz*2)<(this->n*this->n)) ? (fillin+1)*nnz*2 : this->n*this->n;

  ClearMem();

  lindex = 30*this->n+dimfact;          // ?????
  lvalue = dimfact;

  value  = new ARTYPE[lvalue];
  index  = new int[lindex];

} // CreateStructure.
*/

template<class ARTYPE>
inline void ARumSymMatrix<ARTYPE>::ThrowError()
{

  if (status== -1)  {       // Memory is not suficient.
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARumSymMatrix::FactorA");
  }
  else if (status == 1) {    // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARumSymMatrix::FactorA");
  }
  else if (status != 0) {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARumSymMatrix::FactorA");
  }

} // ThrowError.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::FactorA()
{

std::cout <<"ARumSymMatrix::FactorA " << std::endl; 
 
  // Quitting the function if A was not defined.
  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::FactorA");
  }

  ExpandA(); // create Ap Ai Ax

  void *Symbolic ;
  status = umfpack_di_symbolic (this->n, this->n, Ap, Ai, Ax, &Symbolic, NULL, NULL) ;
  ThrowError();
  status =  umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
  ThrowError();
  umfpack_di_free_symbolic (&Symbolic) ;

/*

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to (value, index);

  ExpandA();

  // Decomposing A.

  um2fa(this->n, index[this->n], 0, false, lvalue, lindex, value, 
        index, keep, cntl, icntl, info, rinfo);
*/

  factored = true;

} // FactorA.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::FactorAsI(ARTYPE sigma)
{
std::cout <<"ARumSymMatrix::FactorAsI " << sigma  << std::endl; 

  // Quitting the function if A was not defined.
  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.
  //CreateStructure();

  // Subtracting sigma*I from A.
  ExpandA(sigma);

  // Decomposing AsI.
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  umfpack_di_defaults (Control) ;
  //std::cout << " Ap[n] = " << Ap[this->n] << std::flush;

  void *Symbolic ;
  status = umfpack_di_symbolic (this->n, this->n, Ap, Ai, Ax, &Symbolic, Control, Info) ;
  //std::cout << " symbolic status: " << status << std::endl;
  ThrowError();
  status =  umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL) ;
  //std::cout << " numeric status: " << status << std::endl;
  ThrowError();
  umfpack_di_free_symbolic (&Symbolic) ;

// // Decomposing AsI.
//  um2fa(this->n, index[this->n], 0, false, lvalue, lindex, value,
//        index, keep, cntl, icntl, info, rinfo);


  factored = true;

} // FactorAsI.


template<class ARTYPE>
void ARumSymMatrix<ARTYPE>::MultMv(ARTYPE* v, ARTYPE* w)
{
//std::cout <<"ARumSymMatrix::MultMv ..." << std::flush; 

  int    i,j,k;
  ARTYPE t;

  // Quitting the function if A was not defined.

  if (!this->IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::MultMv");
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
void ARumSymMatrix<ARTYPE>::MultInvv(ARTYPE* v, ARTYPE* w)
{
//std::cout <<"ARumSymMatrix::MultInvv ..." << std::flush; 

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARumSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).

//  ARTYPE* space = new ARTYPE[2*this->n];
//  um2so(this->n, 0, false, lvalue, lindex, value, index,
//        keep, v, w, space, cntl, icntl, info, rinfo);
//  delete[] space;

  status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, w, v, Numeric, NULL, NULL) ;
  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymMatrix::MultInvv");

} // MultInvv.


template<class ARTYPE>
inline void ARumSymMatrix<ARTYPE>::
DefineMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
             int* pcolp, char uplop, double thresholdp,
             int fillinp, bool reducible, bool check)
{

  this->m   = np;
  this->n   = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[this->n]   = nnz;
  uplo      = uplop;
//  fillin    = (fillinp>2) ? fillinp : 2;
  threshold = thresholdp;
//  value     = NULL;
//  index     = NULL;

//  // Preparing umfpack.
//
//  um21i(keep, cntl, icntl, threshold, true, reducible);

  // Checking data.
  if ((check)&&(!DataOK())) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::DefineMatrix");
  }

  this->defined = true;

} // DefineMatrix.


template<class ARTYPE>
inline ARumSymMatrix<ARTYPE>::
ARumSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
              int* pcolp, char uplop, double thresholdp,
              int fillinp, bool reducible, bool check)   : ARMatrix<ARTYPE>(np)
{
  Numeric = NULL;
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;
  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, uplop,
               thresholdp, fillinp, reducible, check);

} // Long constructor.


template<class ARTYPE>
ARumSymMatrix<ARTYPE>::
ARumSymMatrix(const std::string& file, double thresholdp, int fillinp,
              bool reducible, bool check)
{
  Numeric = NULL;
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;

  factored = false;

  try {
    mat.Define(file);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARumSymMatrix");
  }

  if ((mat.NCols() == mat.NRows()) && (mat.IsSymmetric())) {

    DefineMatrix(mat.NCols(), mat.NonZeros(), (ARTYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), 'L', thresholdp,
                 fillinp, reducible, check);
  }
  else {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::ARluSymMatrix");
  }

} // Long constructor (Harwell-Boeing file).


template<class ARTYPE>
ARumSymMatrix<ARTYPE>& ARumSymMatrix<ARTYPE>::
operator=(const ARumSymMatrix<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSMAT_H
