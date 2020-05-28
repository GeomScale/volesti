/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARCSPen.h.
   Arpack++ class ARchSymMPencil definition.
   (CHOLMOD wrapper)

   Author of this class:
      Martin Reuter
      Date 11/05/2012
      
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

#ifndef ARCSPEN_H
#define ARCSPEN_H

//#include "arch.h"
//#include "arerror.h"
#include "blas1c.h"
//#include "lapackc.h"
#include "arcsmat.h"


template<class ARTYPE>
class ARchSymPencil
{

 protected:

  ARchSymMatrix<ARTYPE>* A;
  ARchSymMatrix<ARTYPE>* B;
  cholmod_factor *LAsB ; 
  bool    factoredAsB;
  cholmod_common c ;

  virtual void Copy(const ARchSymPencil& other);

//  void SparseSaxpy(ARTYPE a, ARTYPE x[], int xind[], int nx, ARTYPE y[],
//                   int yind[], int ny, ARTYPE z[], int zind[], int& nz);

//  void ExpandAsB();

//  void SubtractAsB(ARTYPE sigma);

 public:

  bool IsFactored() { return factoredAsB; }

  void FactorAsB(ARTYPE sigma);

  void MultAv(ARTYPE* v, ARTYPE* w) { A->MultMv(v,w); }

  void MultBv(ARTYPE* v, ARTYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(ARTYPE* v, ARTYPE* w);

  void MultInvAsBv(ARTYPE* v, ARTYPE* w);

  void DefineMatrices(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp);

  ARchSymPencil() { factoredAsB = false; A=NULL; B=NULL; LAsB=NULL; cholmod_start (&c) ; }
  // Short constructor that does nothing.

  ARchSymPencil(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp);
  // Long constructor.

  ARchSymPencil(const ARchSymPencil& other) { cholmod_start (&c) ; Copy(other); }
  // Copy constructor.

  virtual ~ARchSymPencil() {  if (LAsB) cholmod_free_factor(&LAsB,&c);  cholmod_finish (&c) ;}
  // Destructor.

  ARchSymPencil& operator=(const ARchSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARchSymPencil member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARTYPE>
inline void ARchSymPencil<ARTYPE>::Copy(const ARchSymPencil<ARTYPE>& other)
{
  if (LAsB) cholmod_free_factor(&LAsB,&c);
  A        = other.A;
  B        = other.B;
  factoredAsB = other.factoredAsB;
  if (factoredAsB)
    LAsB = cholmod_copy_factor(other.LAsB,&c);

} // Copy.

/*
template<class ARTYPE>
void ARchSymPencil<ARTYPE>::
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
void ARchSymPencil<ARTYPE>::ExpandAsB()
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
void ARchSymPencil<ARTYPE>::SubtractAsB(ARTYPE sigma)
{

  int i, acol, bcol, asbcol, scol;

  // Quitting function if A->uplo is not equal to B->uplo.

  if ((A->uplo != B->uplo)&&(sigma != (ARTYPE)0)) {
    throw ArpackError(ArpackError::DIFFERENT_TRIANGLES,
                      "ARchSymPencil::SubtractAsB");
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

} // SubtractAsB.

*/

template<class ARTYPE>
void ARchSymPencil<ARTYPE>::FactorAsB(ARTYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARchSymPencil::FactorAsB");
  }


  if (LAsB) cholmod_free_factor(&LAsB,&c);

  cholmod_sparse* AsB;
  if (sigma != 0.0)
  {
    std::cout << " Subtracting sigma B  (sigma="<<sigma<<")"<<std::endl;
    double alpha[2]; alpha[0]=1.0; alpha[1] = 1.0;
    double beta[2]; beta[0] = -sigma; beta[1]=1.0;
    AsB = cholmod_add(A->A,B->A,alpha,beta,1,0,&c);
  }
  else
    AsB = A->A;
    
//FILE *fp;
//fp=fopen("AsB.asc", "w");
//cholmod_write_sparse(fp,AsB,NULL,NULL,&c);
//FILE *fpa;
//fpa=fopen("As.asc", "w");
//cholmod_write_sparse(fpa,B->A,NULL,NULL,&c);
//FILE *fpb;
//fpb=fopen("Bs.asc", "w");
//cholmod_write_sparse(fpb,A->A,NULL,NULL,&c);

  LAsB = cholmod_analyze (AsB, &c) ;
  int info = cholmod_factorize (AsB, LAsB, &c) ;  

  factoredAsB = (info != 0);  
  if (c.status != CHOLMOD_OK)
  {
    //std::cout << " sigma : " << sigma << std::endl;

    Write_Cholmod_Sparse_Matrix("AsB-error.asc",AsB,&c);

    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARchSymPencil::FactorAsB");
    
    factoredAsB = false;
  }

  if (sigma != 0.0)
    cholmod_free_sparse(&AsB,&c);


} // FactorAsB (ARTYPE shift).


template<class ARTYPE>
void ARchSymPencil<ARTYPE>::MultInvBAv(ARTYPE* v, ARTYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  ::copy(A->ncols(), w, 1, v, 1);
  B->MultInvv(w, w);

} // MultInvBAv.

template<class ARTYPE>
void ARchSymPencil<ARTYPE>::MultInvAsBv(ARTYPE* v, ARTYPE* w)
{
  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARchSymPencil::MultInvAsBv");
  }

  // Solving A.w = v (or AsI.w = v).
  
  //create b from v (data is not copied!!)
  cholmod_dense * b = Create_Cholmod_Dense_Matrix(A->n,1,v,&c);

  cholmod_dense *x = cholmod_solve (CHOLMOD_A, LAsB, b, &c) ;

  Get_Cholmod_Dense_Data(x, A->n, w);

  free(b);
  cholmod_free_dense(&x,&c);


} // MultInvAsBv

template<class ARTYPE>
inline void ARchSymPencil<ARTYPE>::
DefineMatrices(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if (A->n != B->n) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARchSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class ARTYPE>
inline ARchSymPencil<ARTYPE>::
ARchSymPencil(ARchSymMatrix<ARTYPE>& Ap, ARchSymMatrix<ARTYPE>& Bp)
{
  cholmod_start (&c) ;
  LAsB=NULL; 
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class ARTYPE>
ARchSymPencil<ARTYPE>& ARchSymPencil<ARTYPE>::
operator=(const ARchSymPencil<ARTYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSPEN_H
