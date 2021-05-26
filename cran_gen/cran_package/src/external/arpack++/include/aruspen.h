/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARUSPen.h.
   Arpack++ class ARumSymPencil definition.

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

#ifndef ARUSPEN_H
#define ARUSPEN_H

//#include "arch.h"
//#include "arerror.h"
//#include "lapackc.h"
#include "arusmat.h"
#include "blas1c.h"


template<class ARTYPE>
class ARumSymPencil
{

 protected:

  ARumSymMatrix<ARTYPE>* A;
  ARumSymMatrix<ARTYPE>* B;
  //ARumSymMatrix<ARTYPE> AsB;
  void*   Numeric;
  int*    Ap;
  int*    Ai;
  ARTYPE* Ax; 

  virtual void Copy(const ARumSymPencil& other);

//  void SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
//                   int yind[], int ny, ARTYPE z[], int zind[], int& nz);

  void ExpandAsB(ARTYPE sigma);

//  void SubtractAsB(ARTYPE sigma);
  void ClearMem();

 public:

  bool IsFactored() { return (Numeric != NULL); }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  //void MultInvAsBv(ARTYPE* v, ARTYPE* w) { AsB.MultInvv(v,w); }
  void MultInvAsBv(ARTYPE* v, ARTYPE* w);

  void DefineMatrices(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp);

  //ARumSymPencil() { AsB.factored = false; }
  ARumSymPencil() { Numeric = NULL; Ap = NULL; Ai = NULL; Ax = NULL; }
  // Short constructor that does nothing.

  ARumSymPencil(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARumSymPencil(const ARumSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymPencil() { }
  // Destructor.

  ARumSymPencil& operator=(const ARumSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::ClearMem()
{

  if (Numeric) umfpack_di_free_numeric (&Numeric);
  if (Ai) delete [] Ai;
  Ai = NULL;
  if (Ap) delete [] Ap;
  Ap = NULL;
  if (Ax) delete [] Ax;
  Ax = NULL;

} // ClearMem.



template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::Copy(const ARumSymPencil<ARTYPE>& other)
{
  ClearMem();
  A        = other.A;
  B        = other.B;
//  AsB      = other.AsB;

} // Copy.


/*template<class ARTYPE>
void ARumSymPencil<ARTYPE>::
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
void ARumSymPencil<ARTYPE>::ExpandAsB()
{

  int    i, j, k, n;
  int    *pcol, *irow, *index, *pos;
  ARTYPE *value;

  // Initializing variables.

  n     = AsB.n;
  index = AsB.index;
  value = AsB.value;
  irow  = &index[n+1];
  pcol  = new int[AsB.n+1];
  pos   = new int[AsB.n+1];
  for (i=0; i<=n; i++) pcol[i] = index[i];
  for (i=0; i<=n; i++) pos[i] = 0;

  // Counting the elements in each column of AsB.

  if (AsB.uplo == 'U') {

    for (i=0; i!=n; i++) {
      k = pcol[i+1];
      if ((k!=pcol[i])&&(irow[k-1]==i)) k--;
      for (j=pcol[i]; j<k; j++) pos[irow[j]]++;        
    }

  }
  else { // uplo == 'L'

    for (i=0; i!=n; i++) {
      k = pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) k++;
      for (j=k; j<pcol[i+1]; j++) pos[irow[j]]++;        
    }

  }  

  // Summing up index elements.

  for (i=0; i<n; i++) pos[i+1] += pos[i];
  for (i=n; i>0; i--) index[i] += pos[i-1];
    
  // Expanding A.

  if (AsB.uplo == 'U') {

    for (i=n-1; i>=0; i--) {
      pos[i] = index[i]+pcol[i+1]-pcol[i];
      k = pos[i]-1;
      for (j=pcol[i+1]-1; j>=pcol[i]; j--) {
        value[k]  = value[j];
        irow[k--] = irow[j];
      }
    }
    for (i=1; i<n; i++) {
      k = index[i]+pcol[i+1]-pcol[i];
      if ((k>index[i])&&(irow[k-1]==i)) k--;
      for (j=index[i]; j<k; j++) {
        value[pos[irow[j]]]  = value[j];
        irow[pos[irow[j]]++] = i;
      }
    }

  }
  else { // uplo  == 'L'

    for (i=n-1; i>=0; i--) {
      k = index[i+1]-1;
      for (j=pcol[i+1]-1; j>=pcol[i]; j--) {
        value[k]  = value[j];
        irow[k--] = irow[j];
      }
      pos[i] = index[i];
    }
    for (i=0; i<(n-1); i++) {
      k = index[i+1]-pcol[i+1]+pcol[i];
      if ((k<index[i+1])&&(irow[k]==i)) k++;
      for (j=k; j<index[i+1]; j++) {
        value[pos[irow[j]]]  = value[j];
        irow[pos[irow[j]]++] = i;
      }
    }

  }

  AsB.nnz = index[n]; 

  //  Deleting temporary vectors.

  delete[] pcol;  
  delete[] pos;

} // ExpandAsB.


template<class ARTYPE>
void ARumSymPencil<ARTYPE>::SubtractAsB(ARTYPE sigma)
{

  int i, acol, bcol, asbcol, scol;

  // Quitting function if A->uplo is not equal to B->uplo.

  if ((A->uplo != B->uplo)&&(sigma != (ARTYPE)0)) {
    throw ArpackError(ArpackError::DIFFERENT_TRIANGLES,
                      "ARumSymPencil::SubtractAsB");
  }
  AsB.uplo = A->uplo;

  // Subtracting sigma*B from A.

  AsB.index[0] = 0;
  asbcol       = 0;

  for (i=0; i!=AsB.n; i++) {
    bcol = B->pcol[i];
    acol = A->pcol[i];
    SparseSaxpy(-sigma, &B->a[bcol], &B->irow[bcol], B->pcol[i+1]-bcol,
                &A->a[acol], &A->irow[acol], A->pcol[i+1]-acol,
                &AsB.value[asbcol], &AsB.index[asbcol+AsB.n+1], scol);
    asbcol += scol;
    AsB.index[i+1] = asbcol;
  }

  // Expanding AsB.

  ExpandAsB();

  // Adding one to all elements of vector index
  // because the decomposition function was written in FORTRAN.

  for (i=0; i<=AsB.n+AsB.nnz; i++) AsB.index[i]++;

} // SubtractAsB. */


template<class ARTYPE>
void ARumSymPencil<ARTYPE>::ExpandAsB(ARTYPE sigma)
{
std::cout <<"ARumSymPencil::ExpandAsB(" << sigma << ") ..." << std::flush; 

  ClearMem();
 
  int mynnz = 2*A->nnz+2*B->nnz;
  if (sigma == 0.0)
    mynnz = 2*A->nnz;
  
  // create triples (i,j,value)
  int * tripi = new int[mynnz];
  int * tripj = new int[mynnz];
  ARTYPE* tripx = new ARTYPE[mynnz];
  if (tripi == NULL || tripj == NULL || tripx ==NULL)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::ExpandAsB out of memory (1)");
  
  int count = 0;
  int i,j;
  for (i=0; i < A->n; i++)
  {
    // create triplets from A
    for (j=A->pcol[i]; j<(A->pcol[i+1]); j++)
    {
      tripi[count] = i;
      tripj[count] = A->irow[j];
      tripx[count] = A->a[j];
      count++;
      if (i != A->irow[j]) // not on diag
      {
        tripj[count] = i;
        tripi[count] = A->irow[j];
        tripx[count] = A->a[j];
        count++;
      }
    }
  
   if (sigma != 0.0)
   {
    // create triplets from -sigma B
    for (j=B->pcol[i]; j<(B->pcol[i+1]); j++)
    {
      tripi[count] = i;
      tripj[count] = B->irow[j];
      tripx[count] = -sigma * B->a[j];
      count++;
      if (i != B->irow[j]) // not on diag
      {
        tripj[count] = i;
        tripi[count] = B->irow[j];
        tripx[count] = tripx[count-1];
        count++;
      }
    }
    }

  }

  //Write_Triplet_Matrix("A-aruspen.asc",tripi,tripj,tripx,count);

  std::cout<< " ( N = " << A->n << "  NNZ = " << count << " )" << std::flush;
  //std::cout<< " size double " << sizeof(double) << "  size ARTYPE " << sizeof(ARTYPE) << std::endl;
  // convert triples (A-sigma B) to Ax Ap Ai
  Ap = new int[A->n + 1];
  Ai = new int[count];
  Ax = new ARTYPE[count];
  if (!Ap || !Ai || !Ax )
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::ExpandAsB out of memory (2)");
  
  int status = umfpack_di_triplet_to_col (A->n, A->n, count, tripi, tripj, tripx, Ap, Ai, Ax,  (int *)NULL) ;
  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::ExpandAsB triplet to col");

  // cleanup
  delete [] tripi;
  delete [] tripj;
  delete [] tripx;

  //std::cout << std::endl << std::endl;
  //double Control [UMFPACK_CONTROL];
  //Control [UMFPACK_PRL] = 3;
  //status = umfpack_di_report_matrix(A->n, A->n,Ap, Ai, Ax,0,Control);
  //std::cout << " status: " << status << std::endl;
  //std::cout << std::endl << std::endl;

  std::cout <<" done!" << std::endl; 
}

template<class ARTYPE>
void ARumSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumSymPencil::FactorAsB");
  }


  // Subtracting sigma*B from A and storing the result 
  ExpandAsB(sigma);

  // Decomposing AsB.
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  umfpack_di_defaults (Control) ;
  //std::cout <<" loaded defaults" << std::endl;
  void *Symbolic ;
  int status = umfpack_di_symbolic (A->n, A->n, Ap, Ai, Ax, &Symbolic, Control, Info) ;
  std::cout << " symbolic status: " << status << std::endl;
  if (status != UMFPACK_OK)
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB symbolic");
  status =  umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) ;
  std::cout << " numeric status: " << status << std::endl;
  if (status == 1)
  {
    std::cout << " WARNING: MATRIX IS SINGULAR " << std::endl;
    //throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB numeric (matrix singular)");    
  }
  if (status < UMFPACK_OK)
  {
    std::cout << " ERROR CODE: " << status << std::endl;
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB numeric");
  }
  umfpack_di_free_symbolic (&Symbolic) ;

//exit(0);

  // Decomposing AsB.

  //um2fa(AsB.n, AsB.index[AsB.n], 0, false, AsB.lvalue, AsB.lindex, AsB.value,
 //       AsB.index, AsB.keep, AsB.cntl, AsB.icntl, AsB.info, AsB.rinfo);

  // Handling errors.

 // AsB.ThrowError();

 // AsB.factored = true;

} // FactorAsB (ARTYPE shift).


template<class ARTYPE>
void ARumSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.

template<class ARTYPE>
void ARumSymPencil<ARTYPE>::MultInvAsBv(ARTYPE* v, ARTYPE* w)
{
  if (!Numeric) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARchSymPencil::MultInvAsBv");
  }

  // Solving A.w = v (or AsI.w = v).
   int status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, w, v, Numeric, NULL, NULL) ;
  if (status == 1)
  {
    std::cout << " WARNING: MATRIX IS SINGULAR " << std::endl;
    //throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::FactorAsB numeric (matrix singular)");    
  }
  if (status < UMFPACK_OK)
  {
    std::cout << " ERROR CODE: " << status << std::endl;
    throw ArpackError(ArpackError::PARAMETER_ERROR, "ARumSymPencil::MultInvAsBv");
 
  }

} // MultInvAsBv

template<class ARTYPE>
inline void ARumSymPencil<ARTYPE>::
DefineMatrices(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if (A->n != B->n) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARumSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE>
inline ARumSymPencil<ARTYPE>::
ARumSymPencil(ARumSymMatrix<ARTYPE>& Ap, ARumSymMatrix<ARTYPE>& Bp)
{
  Numeric = NULL;
  Ap = NULL;
  Ai = NULL;
  Ax = NULL;

  //AsB.factored  = false;
  DefineMatrices(Ap, Bp);
  

} // Long constructor.


template<class ARTYPE>
ARumSymPencil<ARTYPE>& ARumSymPencil<ARTYPE>::
operator=(const ARumSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSPEN_H
