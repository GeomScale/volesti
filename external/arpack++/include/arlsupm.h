/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARLSupM.h.
   Unaltered copy of supermatrix.h (from SuperLU package).
*/

#ifndef __SUPERLU_SUPERMATRIX /* allow multiple inclusions */
#define __SUPERLU_SUPERMATRIX


/********************************************
 * The matrix types are defined as follows. *
 ********************************************/
typedef enum {
    SLU_NC,    /* column-wise, no supernode */
    SLU_NCP,   /* column-wise, column-permuted, no supernode 
                  (The consecutive columns of nonzeros, after permutation,
		   may not be stored  contiguously.) */
    SLU_NR,    /* row-wize, no supernode */
    SLU_SC,    /* column-wise, supernode */
    SLU_SCP,   /* supernode, column-wise, permuted */    
    SLU_SR,    /* row-wise, supernode */
    SLU_DN,     /* Fortran style column-wise storage for dense matrix */
    SLU_NR_loc  /* distributed compressed row format  */ 
} Stype_t;

typedef enum {
    SLU_S,     /* single */
    SLU_D,     /* double */
    SLU_C,     /* single complex */
    SLU_Z      /* double complex */
} Dtype_t;

typedef enum {
    SLU_GE,    /* general */
    SLU_TRLU,  /* lower triangular, unit diagonal */
    SLU_TRUU,  /* upper triangular, unit diagonal */
    SLU_TRL,   /* lower triangular */
    SLU_TRU,   /* upper triangular */
    SLU_SYL,   /* symmetric, store lower half */
    SLU_SYU,   /* symmetric, store upper half */
    SLU_HEL,   /* Hermitian, store lower half */
    SLU_HEU    /* Hermitian, store upper half */
} Mtype_t;

typedef struct {
	Stype_t Stype; /* Storage type: interprets the storage structure 
		   	  pointed to by *Store. */
	Dtype_t Dtype; /* Data type. */
	Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
			  the matrix. */
	int  nrow;   /* number of rows */
	int  ncol;   /* number of columns */
	void *Store;   /* pointer to the actual storage of the matrix */
} SuperMatrix;

/***********************************************
 * The storage schemes are defined as follows. *
 ***********************************************/

/* Stype == SLU_NC (Also known as Harwell-Boeing sparse matrix format) */
typedef struct {
    int  nnz;	    /* number of nonzeros in the matrix */
    void *nzval;    /* pointer to array of nonzero values, packed by column */
    int  *rowind; /* pointer to array of row indices of the nonzeros */
    int  *colptr; /* pointer to array of beginning of columns in nzval[] 
		       and rowind[]  */
                    /* Note:
		       Zero-based indexing is used;
		       colptr[] has ncol+1 entries, the last one pointing
		       beyond the last column, so that colptr[ncol] = nnz. */
} NCformat;

/* Stype == SLU_NR */
typedef struct {
    int  nnz;	    /* number of nonzeros in the matrix */
    void *nzval;    /* pointer to array of nonzero values, packed by raw */
    int  *colind; /* pointer to array of columns indices of the nonzeros */
    int  *rowptr; /* pointer to array of beginning of rows in nzval[] 
		       and colind[]  */
                    /* Note:
		       Zero-based indexing is used;
		       rowptr[] has nrow+1 entries, the last one pointing
		       beyond the last row, so that rowptr[nrow] = nnz. */
} NRformat;

/* Stype == SLU_SC */
typedef struct {
  int  nnz;	     /* number of nonzeros in the matrix */
  int  nsuper;     /* number of supernodes, minus 1 */
  void *nzval;       /* pointer to array of nonzero values, packed by column */
  int *nzval_colptr;/* pointer to array of beginning of columns in nzval[] */
  int *rowind;     /* pointer to array of compressed row indices of 
			rectangular supernodes */
  int *rowind_colptr;/* pointer to array of beginning of columns in rowind[] */
  int *col_to_sup;   /* col_to_sup[j] is the supernode number to which column 
			j belongs; mapping from column to supernode number. */
  int *sup_to_col;   /* sup_to_col[s] points to the start of the s-th 
			supernode; mapping from supernode number to column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
		              sup_to_col: 0 1 2 4 7 12           (nsuper=4) */
                     /* Note:
		        Zero-based indexing is used;
		        nzval_colptr[], rowind_colptr[], col_to_sup and
		        sup_to_col[] have ncol+1 entries, the last one
		        pointing beyond the last column.
		        For col_to_sup[], only the first ncol entries are
		        defined. For sup_to_col[], only the first nsuper+2
		        entries are defined. */
} SCformat;

/* Stype == SLU_SCP */
typedef struct {
  int  nnz;	     /* number of nonzeros in the matrix */
  int  nsuper;     /* number of supernodes */
  void *nzval;       /* pointer to array of nonzero values, packed by column */
  int  *nzval_colbeg;/* nzval_colbeg[j] points to beginning of column j
			  in nzval[] */
  int  *nzval_colend;/* nzval_colend[j] points to one past the last element
			  of column j in nzval[] */
  int  *rowind;      /* pointer to array of compressed row indices of 
			  rectangular supernodes */
  int *rowind_colbeg;/* rowind_colbeg[j] points to beginning of column j
			  in rowind[] */
  int *rowind_colend;/* rowind_colend[j] points to one past the last element
			  of column j in rowind[] */
  int *col_to_sup;   /* col_to_sup[j] is the supernode number to which column
			  j belongs; mapping from column to supernode. */
  int *sup_to_colbeg; /* sup_to_colbeg[s] points to the start of the s-th 
			   supernode; mapping from supernode to column.*/
  int *sup_to_colend; /* sup_to_colend[s] points to one past the end of the
			   s-th supernode; mapping from supernode number to
			   column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
		              sup_to_colbeg: 0 1 2 4 7              (nsuper=4)
			      sup_to_colend: 1 2 4 7 12                    */
                     /* Note:
		        Zero-based indexing is used;
		        nzval_colptr[], rowind_colptr[], col_to_sup and
		        sup_to_col[] have ncol+1 entries, the last one
		        pointing beyond the last column.         */
} SCPformat;

/* Stype == SLU_NCP */
typedef struct {
    int nnz;	  /* number of nonzeros in the matrix */
    void *nzval;  /* pointer to array of nonzero values, packed by column */
    int *rowind;/* pointer to array of row indices of the nonzeros */
		  /* Note: nzval[]/rowind[] always have the same length */
    int *colbeg;/* colbeg[j] points to the beginning of column j in nzval[] 
                     and rowind[]  */
    int *colend;/* colend[j] points to one past the last element of column
		     j in nzval[] and rowind[]  */
		  /* Note:
		     Zero-based indexing is used;
		     The consecutive columns of the nonzeros may not be 
		     contiguous in storage, because the matrix has been 
		     postmultiplied by a column permutation matrix. */
} NCPformat;

/* Stype == SLU_DN */
typedef struct {
    int lda;    /* leading dimension */
    void *nzval;  /* array of size lda*ncol to represent a dense matrix */
} DNformat;

/* Stype == SLU_NR_loc (Distributed Compressed Row Format) */
typedef struct {
    int nnz_loc;   /* number of nonzeros in the local submatrix */
    int m_loc;     /* number of rows local to this processor */
    int fst_row;   /* global index of the first row */
    void  *nzval;    /* pointer to array of nonzero values, packed by row */
    int *rowptr;   /* pointer to array of beginning of rows in nzval[] 
			and colind[]  */
    int *colind;   /* pointer to array of column indices of the nonzeros */
                     /* Note:
			Zero-based indexing is used;
			rowptr[] has n_loc + 1 entries, the last one pointing
			beyond the last row, so that rowptr[n_loc] = nnz_loc.*/
} NRformat_loc;


#endif  /* __SUPERLU_SUPERMATRIX */
