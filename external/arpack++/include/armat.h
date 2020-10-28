/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARMat.h
   Generic matrix template with a matrix-vector product.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARMAT_H
#define ARMAT_H

template<class ARTYPE>
class ARMatrix {

 protected:

  int  m, n;    // Number of rows and columns.
  bool defined;
 
 public:

  ARMatrix() { defined = false; }
  // Short constructor.

  ARMatrix(int nrows, int ncols = 0)
  // Long constructor.
  {
    m = nrows;
    n = (ncols?ncols:nrows);
    defined = false;
  } // Constructor.

  virtual ~ARMatrix() { }
  // Destructor.

  int nrows() { return m; }

  int ncols() { return n; }

  bool IsDefined() { return defined; }

  virtual void MultMv(ARTYPE* v, ARTYPE* w) = 0;
  // Matrix-vector product: w = A*v.

}; // ARMatrix.

#endif // ARMAT_H

