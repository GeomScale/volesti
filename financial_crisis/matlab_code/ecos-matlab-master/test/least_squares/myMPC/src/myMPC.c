/*
FORCES - Fast interior point code generation for multistage problems.
Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
Automatic Control Laboratory, ETH Zurich.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../include/myMPC.h"

/* for square root */
#include <math.h> 

/* SYSTEM INCLUDES FOR PRINTING ---------------------------------------- */
#ifndef USEMEXPRINTS
#include <stdio.h>
#define PRINTTEXT printf
#else
#include "mex.h"
#define PRINTTEXT mexPrintf
#endif

/* TIMING LIBRARY ------------------------------------------------- */

/* ARE WE ON WINDOWS? */
#if (defined WIN32 || defined _WIN64 || defined _WIN32)

/* Use Windows QueryPerformanceCounter for timing */

#include <Windows.h>

typedef struct myMPC_timer{
	LARGE_INTEGER tic;
	LARGE_INTEGER toc;
	LARGE_INTEGER freq;
} myMPC_timer;


void myMPC_tic(myMPC_timer* t)
{
	QueryPerformanceFrequency(&t->freq);
	QueryPerformanceCounter(&t->tic);
}



myMPC_FLOAT myMPC_toc(myMPC_timer* t)
{
	QueryPerformanceCounter(&t->toc);
	return ((t->toc.QuadPart - t->tic.QuadPart) / (myMPC_FLOAT)t->freq.QuadPart);
}


/* WE ARE ON THE MAC */
#elif (defined __APPLE__)
#include <mach/mach_time.h>


/* Use MAC OSX  mach_time for timing */
typedef struct myMPC_timer{
	uint64_t tic;
	uint64_t toc;
	mach_timebase_info_data_t tinfo;

} myMPC_timer;


void myMPC_tic(myMPC_timer* t)
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}



myMPC_FLOAT myMPC_toc(myMPC_timer* t)
{
    uint64_t duration; /* elapsed time in clock cycles*/
    t->toc = mach_absolute_time();
	duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (myMPC_FLOAT)duration / 1000000000;
}



/* WE ARE ON SOME OTHER UNIX/LINUX SYSTEM */
#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <time.h>
typedef struct myMPC_timer{
	struct timespec tic;
	struct timespec toc;
} myMPC_timer;


/* read current time */
void myMPC_tic(myMPC_timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}



/* return time passed since last call to tic on this timer */
double myMPC_toc(myMPC_timer* t)
{
	struct timespec temp;
	clock_gettime(CLOCK_MONOTONIC, &t->toc);	

	if ((t->toc.tv_nsec - t->tic.tv_nsec)<0) {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec-1;
		temp.tv_nsec = 1000000000+t->toc.tv_nsec - t->tic.tv_nsec;
	} else {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}

	return (myMPC_FLOAT)temp.tv_sec + (myMPC_FLOAT)temp.tv_nsec / 1000000000;
}


#endif

/* LINEAR ALGEBRA LIBRARY ---------------------------------------------- */
/*
 * Initializes a vector of length 29 with a value.
 */
void myMPC_LA_INITIALIZEVECTOR_29(myMPC_FLOAT* vec, myMPC_FLOAT value)
{
	int i;
	for( i=0; i<29; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 20 with a value.
 */
void myMPC_LA_INITIALIZEVECTOR_20(myMPC_FLOAT* vec, myMPC_FLOAT value)
{
	int i;
	for( i=0; i<20; i++ )
	{
		vec[i] = value;
	}
}


/*
 * Initializes a vector of length 58 with a value.
 */
void myMPC_LA_INITIALIZEVECTOR_58(myMPC_FLOAT* vec, myMPC_FLOAT value)
{
	int i;
	for( i=0; i<58; i++ )
	{
		vec[i] = value;
	}
}


/* 
 * Calculates a dot product and adds it to a variable: z += x'*y; 
 * This function is for vectors of length 58.
 */
void myMPC_LA_DOTACC_58(myMPC_FLOAT *x, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<58; i++ ){
		*z += x[i]*y[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, dense matrix of size [3 x 3]
 *             f  - column vector of size 3
 *             z  - column vector of size 3
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 3
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void myMPC_LA_DENSE_QUADFCN_3(myMPC_FLOAT* H, myMPC_FLOAT* f, myMPC_FLOAT* z, myMPC_FLOAT* grad, myMPC_FLOAT* value)
{
	int i;
	int j;
	int k = 0;
	myMPC_FLOAT hz;	
	for( i=0; i<3; i++){
		hz = 0;
		for( j=0; j<3; j++ )
		{
			hz += H[k++]*z[j];
		}
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/*
 * Calculates the gradient and the value for a quadratic function 0.5*z'*H*z + f'*z
 *
 * INPUTS:     H  - Symmetric Hessian, dense matrix of size [2 x 2]
 *             f  - column vector of size 2
 *             z  - column vector of size 2
 *
 * OUTPUTS: grad  - gradient at z (= H*z + f), column vector of size 2
 *          value <-- value + 0.5*z'*H*z + f'*z (value will be modified)
 */
void myMPC_LA_DENSE_QUADFCN_2(myMPC_FLOAT* H, myMPC_FLOAT* f, myMPC_FLOAT* z, myMPC_FLOAT* grad, myMPC_FLOAT* value)
{
	int i;
	int j;
	int k = 0;
	myMPC_FLOAT hz;	
	for( i=0; i<2; i++){
		hz = 0;
		for( j=0; j<2; j++ )
		{
			hz += H[k++]*z[j];
		}
		grad[i] = hz + f[i];
		*value += 0.5*hz*z[i] + f[i]*z[i];
	}
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void myMPC_LA_DENSE_2MVMSUB_4_3_3(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *b, myMPC_FLOAT *l, myMPC_FLOAT *r, myMPC_FLOAT *z, myMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	myMPC_FLOAT AxBu[4];
	myMPC_FLOAT norm = *y;
	myMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<4; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<4; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<3; n++ ){
		for( i=0; i<4; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<4; i++ ){
		r[i] = AxBu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *b, myMPC_FLOAT *l, myMPC_FLOAT *r, myMPC_FLOAT *z, myMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	myMPC_FLOAT AxBu[2];
	myMPC_FLOAT norm = *y;
	myMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<2; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<3; n++ ){
		for( i=0; i<2; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<2; i++ ){
		r[i] = AxBu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/* 
 * Computes r = A*x + B*u - b
 * and      y = max([norm(r,inf), y])
 * and      z -= l'*r
 * where A is stored in column major format
 */
void myMPC_LA_DENSE_2MVMSUB_2_3_2(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *b, myMPC_FLOAT *l, myMPC_FLOAT *r, myMPC_FLOAT *z, myMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;
	myMPC_FLOAT AxBu[2];
	myMPC_FLOAT norm = *y;
	myMPC_FLOAT lr = 0;

	/* do A*x + B*u first */
	for( i=0; i<2; i++ ){
		AxBu[i] = A[k++]*x[0] + B[m++]*u[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<2; i++ ){
			AxBu[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<2; n++ ){
		for( i=0; i<2; i++ ){
			AxBu[i] += B[m++]*u[n];
		}		
	}

	for( i=0; i<2; i++ ){
		r[i] = AxBu[i] - b[i];
		lr += l[i]*r[i];
		if( r[i] > norm ){
			norm = r[i];
		}
		if( -r[i] > norm ){
			norm = -r[i];
		}
	}
	*y = norm;
	*z -= lr;
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [4 x 3]
 * and stored in column major format. Note the transpose of M!
 */
void myMPC_LA_DENSE_MTVM_4_3(myMPC_FLOAT *M, myMPC_FLOAT *x, myMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<3; i++ ){
		y[i] = 0;
		for( j=0; j<4; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [2 x 3]
 * and B is of size [4 x 3]
 * and stored in column major format. Note the transposes of A and B!
 */
void myMPC_LA_DENSE_MTVM2_2_3_4(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<3; i++ ){
		z[i] = 0;
		for( j=0; j<2; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<4; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication z = A'*x + B'*y 
 * where A is of size [2 x 3]
 * and B is of size [2 x 3]
 * and stored in column major format. Note the transposes of A and B!
 */
void myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	int j;
	int k = 0;
	int n;
	int m = 0;
	for( i=0; i<3; i++ ){
		z[i] = 0;
		for( j=0; j<2; j++ ){
			z[i] += A[k++]*x[j];
		}
		for( n=0; n<2; n++ ){
			z[i] += B[m++]*y[n];
		}
	}
}


/*
 * Matrix vector multiplication y = M'*x where M is of size [2 x 2]
 * and stored in column major format. Note the transpose of M!
 */
void myMPC_LA_DENSE_MTVM_2_2(myMPC_FLOAT *M, myMPC_FLOAT *x, myMPC_FLOAT *y)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<2; i++ ){
		y[i] = 0;
		for( j=0; j<2; j++ ){
			y[i] += M[k++]*x[j];
		}
	}
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 3. Output z is of course scalar.
 */
void myMPC_LA_VSUBADD3_3(myMPC_FLOAT* t, myMPC_FLOAT* u, int* uidx, myMPC_FLOAT* v, myMPC_FLOAT* w, myMPC_FLOAT* y, myMPC_FLOAT* z, myMPC_FLOAT* r)
{
	int i;
	myMPC_FLOAT norm = *r;
	myMPC_FLOAT vx = 0;
	myMPC_FLOAT x;
	for( i=0; i<3; i++){
		x = t[i] - u[uidx[i]];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 3. Output z is of course scalar.
 */
void myMPC_LA_VSUBADD2_3(myMPC_FLOAT* t, int* tidx, myMPC_FLOAT* u, myMPC_FLOAT* v, myMPC_FLOAT* w, myMPC_FLOAT* y, myMPC_FLOAT* z, myMPC_FLOAT* r)
{
	int i;
	myMPC_FLOAT norm = *r;
	myMPC_FLOAT vx = 0;
	myMPC_FLOAT x;
	for( i=0; i<3; i++){
		x = t[tidx[i]] - u[i];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 2. Output z is of course scalar.
 */
void myMPC_LA_VSUBADD3_2(myMPC_FLOAT* t, myMPC_FLOAT* u, int* uidx, myMPC_FLOAT* v, myMPC_FLOAT* w, myMPC_FLOAT* y, myMPC_FLOAT* z, myMPC_FLOAT* r)
{
	int i;
	myMPC_FLOAT norm = *r;
	myMPC_FLOAT vx = 0;
	myMPC_FLOAT x;
	for( i=0; i<2; i++){
		x = t[i] - u[uidx[i]];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Vector subtraction and addition.
 *	 Input: five vectors t, tidx, u, v, w and two scalars z and r
 *	 Output: y = t(tidx) - u + w
 *           z = z - v'*x;
 *           r = max([norm(y,inf), z]);
 * for vectors of length 2. Output z is of course scalar.
 */
void myMPC_LA_VSUBADD2_2(myMPC_FLOAT* t, int* tidx, myMPC_FLOAT* u, myMPC_FLOAT* v, myMPC_FLOAT* w, myMPC_FLOAT* y, myMPC_FLOAT* z, myMPC_FLOAT* r)
{
	int i;
	myMPC_FLOAT norm = *r;
	myMPC_FLOAT vx = 0;
	myMPC_FLOAT x;
	for( i=0; i<2; i++){
		x = t[tidx[i]] - u[i];
		y[i] = x + w[i];
		vx += v[i]*x;
		if( y[i] > norm ){
			norm = y[i];
		}
		if( -y[i] > norm ){
			norm = -y[i];
		}
	}
	*z -= vx;
	*r = norm;
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 3
 * Returns also L/S, a value that is often used elsewhere.
 */
void myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_FLOAT *lu, myMPC_FLOAT *su, myMPC_FLOAT *ru, myMPC_FLOAT *ll, myMPC_FLOAT *sl, myMPC_FLOAT *rl, int* lbIdx, int* ubIdx, myMPC_FLOAT *grad, myMPC_FLOAT *lubysu, myMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<3; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<3; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<3; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Computes inequality constraints gradient-
 * Special function for box constraints of length 2
 * Returns also L/S, a value that is often used elsewhere.
 */
void myMPC_LA_INEQ_B_GRAD_2_2_2(myMPC_FLOAT *lu, myMPC_FLOAT *su, myMPC_FLOAT *ru, myMPC_FLOAT *ll, myMPC_FLOAT *sl, myMPC_FLOAT *rl, int* lbIdx, int* ubIdx, myMPC_FLOAT *grad, myMPC_FLOAT *lubysu, myMPC_FLOAT *llbysl)
{
	int i;
	for( i=0; i<2; i++ ){
		grad[i] = 0;
	}
	for( i=0; i<2; i++ ){		
		llbysl[i] = ll[i] / sl[i];
		grad[lbIdx[i]] -= llbysl[i]*rl[i];
	}
	for( i=0; i<2; i++ ){
		lubysu[i] = lu[i] / su[i];
		grad[ubIdx[i]] += lubysu[i]*ru[i];
	}
}


/*
 * Addition of three vectors  z = u + w + v
 * of length 29.
 */
void myMPC_LA_VVADD3_29(myMPC_FLOAT *u, myMPC_FLOAT *v, myMPC_FLOAT *w, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<29; i++ ){
		z[i] = u[i] + v[i] + w[i];
	}
}


/*
 * Special function to compute the Dense positive definite 
 * augmented Hessian for block size 3.
 *
 * Inputs: - H = dense cost Hessian in column major storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = H + diag(llbysl) + diag(lubysu)
 * where Phi is stored in lower triangular row major format
 */
void myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_FLOAT *H, myMPC_FLOAT *llbysl, int* lbIdx, myMPC_FLOAT *lubysu, int* ubIdx, myMPC_FLOAT *Phi)
{
	int i;
	int j;
	int k = 0;
	
	/* copy lower triangular part of H into PHI */
	for( i=0; i<3; i++ ){
		for( j=0; j<=i; j++ ){
			Phi[k++] = H[i*3+j];
		}		
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<3; i++ ){
		j = lbIdx[i];
		Phi[((j+1)*(j+2))/2-1] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<3; i++){
		j = ubIdx[i];
		Phi[((j+1)*(j+2))/2-1] +=  lubysu[i];
	}

}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 3.
 */
void myMPC_LA_DENSE_CHOL2_3(myMPC_FLOAT *A)
{
    int i, j, k, ii, jj, di, dj;
    myMPC_FLOAT l;
    myMPC_FLOAT Mii;
    
	ii=0; di=0;
    for( i=0; i<3; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += A[ii+k]*A[ii+k];
        }        
        
        Mii = A[ii+i] - l;
        
        if( Mii < 1e-13 ){
#if myMPC_SET_PRINTLEVEL > 0
#ifdef PRINTNUMERICALWARNINGS
            PRINTTEXT("WARNING: pivot in Cholesky factorization close to 0, regularizing...\n");
#endif
#endif
            Mii = 4e-4;
        }
            
        A[ii+i] = sqrt(Mii);        

		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<3; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += A[jj+k]*A[ii+k];
            }
            A[jj+i] = (A[jj+i] - l)/A[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [4 x 3],
 * B is given and of size [4 x 3], L is a lower tri-
 * angular matrix of size 3 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void myMPC_LA_DENSE_MATRIXFORWARDSUB_4_3(myMPC_FLOAT *L, myMPC_FLOAT *B, myMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    myMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<3; j++ ){        
        for( i=0; i<4; i++ ){
            a = B[j*4+i];
            for( k=0; k<j; k++ ){
                a -= A[k*4+i]*L[ii+k];
            }
            A[j*4+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [4 x 3] in column
 * storage format, and B is of size [4 x 3] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void myMPC_LA_DENSE_MMT2_4_3_3(myMPC_FLOAT *A, myMPC_FLOAT *B, myMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    myMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*4+i]*A[k*4+j];
            }			
			for( k=0; k<3; k++ ){
                ltemp += B[k*4+i]*B[k*4+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 4 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void myMPC_LA_DENSE_CHOL_4(myMPC_FLOAT *A, myMPC_FLOAT *L)
{
    int i, j, k, ii, jj, di, dj;
    myMPC_FLOAT l;
    myMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<4; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<4; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
        if( Mii < 1e-13 ){
#if myMPC_SET_PRINTLEVEL > 0
#ifdef PRINTNUMERICALWARNINGS
            PRINTTEXT("WARNING: pivot in Cholesky factorization close to 0, regularizing...\n");
#endif
#endif
            Mii = 4e-4;
        }
            
        L[ii+i] = sqrt(Mii);        

		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<4; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += L[jj+k]*L[ii+k];
            }
            L[jj+i] = (L[jj+i] - l)/L[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }	
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 3.
 */
void myMPC_LA_DENSE_FORWARDSUB_3(myMPC_FLOAT *L, myMPC_FLOAT *b, myMPC_FLOAT *y)
{
    int i,j,ii,di;
    myMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<3; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }
        y[i] = yel / L[ii+i];
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A an B are stored in column major format
 */
void myMPC_LA_DENSE_2MVMSUB2_4_3_3(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<4; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<3; n++ ){
		for( i=0; i<4; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 4.
 */
void myMPC_LA_DENSE_FORWARDSUB_4(myMPC_FLOAT *L, myMPC_FLOAT *b, myMPC_FLOAT *y)
{
    int i,j,ii,di;
    myMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<4; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }
        y[i] = yel / L[ii+i];
        ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [2 x 3],
 * B is given and of size [2 x 3], L is a lower tri-
 * angular matrix of size 3 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_FLOAT *L, myMPC_FLOAT *B, myMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    myMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<3; j++ ){        
        for( i=0; i<2; i++ ){
            a = B[j*2+i];
            for( k=0; k<j; k++ ){
                a -= A[k*2+i]*L[ii+k];
            }
            A[j*2+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [2 x 3] in column
 * storage format, and B is of size [2 x 3] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void myMPC_LA_DENSE_MMT2_2_3_3(myMPC_FLOAT *A, myMPC_FLOAT *B, myMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    myMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
            }			
			for( k=0; k<3; k++ ){
                ltemp += B[k*2+i]*B[k*2+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [4 x 3]
 *  size(B) = [2 x 3]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void myMPC_LA_DENSE_MMTM_4_3_2(myMPC_FLOAT *A, myMPC_FLOAT *B, myMPC_FLOAT *C)
{
    int i, j, k;
    myMPC_FLOAT temp;
    
    for( i=0; i<4; i++ ){        
        for( j=0; j<2; j++ ){
            temp = 0; 
            for( k=0; k<3; k++ ){
                temp += A[k*4+i]*B[k*2+j];
            }						
            C[j*4+i] = temp;
        }
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B'
 * where A is to be computed and is of size [2 x 4],
 * B is given and of size [2 x 4], L is a lower tri-
 * angular matrix of size 4 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_4(myMPC_FLOAT *L, myMPC_FLOAT *B, myMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    myMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<4; j++ ){        
        for( i=0; i<2; i++ ){
            a = B[i*4+j];
            for( k=0; k<j; k++ ){
                a -= A[k*2+i]*L[ii+k];
            }            
			A[j*2+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 2
 * and A is a dense matrix of size [2 x 4] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void myMPC_LA_DENSE_MMTSUB_2_4(myMPC_FLOAT *A, myMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    myMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<4; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 2 and outputting
 * the Cholesky factor to matrix L in lower triangular format.
 */
void myMPC_LA_DENSE_CHOL_2(myMPC_FLOAT *A, myMPC_FLOAT *L)
{
    int i, j, k, ii, jj, di, dj;
    myMPC_FLOAT l;
    myMPC_FLOAT Mii;

	/* copy A to L first and then operate on L */
	/* COULD BE OPTIMIZED */
	ii=0; di=0;
	for( i=0; i<2; i++ ){
		for( j=0; j<=i; j++ ){
			L[ii+j] = A[ii+j];
		}
		ii += ++di;
	}    
	
	/* factor L */
	ii=0; di=0;
    for( i=0; i<2; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += L[ii+k]*L[ii+k];
        }        
        
        Mii = L[ii+i] - l;
        
        if( Mii < 1e-13 ){
#if myMPC_SET_PRINTLEVEL > 0
#ifdef PRINTNUMERICALWARNINGS
            PRINTTEXT("WARNING: pivot in Cholesky factorization close to 0, regularizing...\n");
#endif
#endif
            Mii = 4e-4;
        }
            
        L[ii+i] = sqrt(Mii);        

		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<2; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += L[jj+k]*L[ii+k];
            }
            L[jj+i] = (L[jj+i] - l)/L[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }	
}


/* 
 * Computes r = b - A*x - B*u
 * where A an B are stored in column major format
 */
void myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<3; n++ ){
		for( i=0; i<2; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = b - A*x
 * where A is stored in column major format
 */
void myMPC_LA_DENSE_MVMSUB1_2_4(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<4; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Forward substitution to solve L*y = b where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * The dimensions involved are 2.
 */
void myMPC_LA_DENSE_FORWARDSUB_2(myMPC_FLOAT *L, myMPC_FLOAT *b, myMPC_FLOAT *y)
{
    int i,j,ii,di;
    myMPC_FLOAT yel;
            
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }
        y[i] = yel / L[ii+i];
        ii += ++di;
    }
}


/**
 * Compute C = A*B' where 
 *
 *	size(A) = [2 x 3]
 *  size(B) = [2 x 3]
 * 
 * and all matrices are stored in column major format.
 *
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE.  
 * 
 */
void myMPC_LA_DENSE_MMTM_2_3_2(myMPC_FLOAT *A, myMPC_FLOAT *B, myMPC_FLOAT *C)
{
    int i, j, k;
    myMPC_FLOAT temp;
    
    for( i=0; i<2; i++ ){        
        for( j=0; j<2; j++ ){
            temp = 0; 
            for( k=0; k<3; k++ ){
                temp += A[k*2+i]*B[k*2+j];
            }						
            C[j*2+i] = temp;
        }
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B'
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a lower tri-
 * angular matrix of size 2 stored in lower triangular 
 * storage format. Note the transpose of L AND B!
 *
 * Result: A in column major storage format.
 *
 */
void myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_FLOAT *L, myMPC_FLOAT *B, myMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    myMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<2; j++ ){        
        for( i=0; i<2; i++ ){
            a = B[i*2+j];
            for( k=0; k<j; k++ ){
                a -= A[k*2+i]*L[ii+k];
            }            
			A[j*2+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Compute L = L - A*A', where L is lower triangular of size 2
 * and A is a dense matrix of size [2 x 2] in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void myMPC_LA_DENSE_MMTSUB_2_2(myMPC_FLOAT *A, myMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    myMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<2; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
            }						
            L[ii+j] -= ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x
 * where A is stored in column major format
 */
void myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0];
	}	
	for( j=1; j<2; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Special function to compute the Dense positive definite 
 * augmented Hessian for block size 2.
 *
 * Inputs: - H = dense cost Hessian in column major storage format
 *         - llbysl = L / S of lower bounds
 *         - lubysu = L / S of upper bounds
 *
 * Output: Phi = H + diag(llbysl) + diag(lubysu)
 * where Phi is stored in lower triangular row major format
 */
void myMPC_LA_INEQ_DENSE_HESS_2_2_2(myMPC_FLOAT *H, myMPC_FLOAT *llbysl, int* lbIdx, myMPC_FLOAT *lubysu, int* ubIdx, myMPC_FLOAT *Phi)
{
	int i;
	int j;
	int k = 0;
	
	/* copy lower triangular part of H into PHI */
	for( i=0; i<2; i++ ){
		for( j=0; j<=i; j++ ){
			Phi[k++] = H[i*2+j];
		}		
	}

	/* add llbysl onto Phi where necessary */
	for( i=0; i<2; i++ ){
		j = lbIdx[i];
		Phi[((j+1)*(j+2))/2-1] += llbysl[i];
	}

	/* add lubysu onto Phi where necessary */
	for( i=0; i<2; i++){
		j = ubIdx[i];
		Phi[((j+1)*(j+2))/2-1] +=  lubysu[i];
	}

}


/**
 * Cholesky factorization as above, but working on a matrix in 
 * lower triangular storage format of size 2.
 */
void myMPC_LA_DENSE_CHOL2_2(myMPC_FLOAT *A)
{
    int i, j, k, ii, jj, di, dj;
    myMPC_FLOAT l;
    myMPC_FLOAT Mii;
    
	ii=0; di=0;
    for( i=0; i<2; i++ ){
        l = 0;
        for( k=0; k<i; k++ ){
            l += A[ii+k]*A[ii+k];
        }        
        
        Mii = A[ii+i] - l;
        
        if( Mii < 1e-13 ){
#if myMPC_SET_PRINTLEVEL > 0
#ifdef PRINTNUMERICALWARNINGS
            PRINTTEXT("WARNING: pivot in Cholesky factorization close to 0, regularizing...\n");
#endif
#endif
            Mii = 4e-4;
        }
            
        A[ii+i] = sqrt(Mii);        

		jj = ((i+1)*(i+2))/2; dj = i+1;
        for( j=i+1; j<2; j++ ){
            l = 0;            
            for( k=0; k<i; k++ ){
                l += A[jj+k]*A[ii+k];
            }
            A[jj+i] = (A[jj+i] - l)/A[ii+i];            
			jj += ++dj;
        }
		ii += ++di;
    }
}


/**
 * Forward substitution for the matrix equation A*L' = B
 * where A is to be computed and is of size [2 x 2],
 * B is given and of size [2 x 2], L is a lower tri-
 * angular matrix of size 2 stored in lower triangular 
 * storage format. Note the transpose of L!
 *
 * Result: A in column major storage format.
 *
 */
void myMPC_LA_DENSE_MATRIXFORWARDSUB_2_2(myMPC_FLOAT *L, myMPC_FLOAT *B, myMPC_FLOAT *A)
{
    int i,j,k,ii,di;
    myMPC_FLOAT a;
    
    ii=0; di=0;
    for( j=0; j<2; j++ ){        
        for( i=0; i<2; i++ ){
            a = B[j*2+i];
            for( k=0; k<j; k++ ){
                a -= A[k*2+i]*L[ii+k];
            }
            A[j*2+i] = a/L[ii+j];
        }
        ii += ++di;
    }
}


/**
 * Compute L = A*A' + B*B', where L is lower triangular of size NXp1
 * and A is a dense matrix of size [2 x 3] in column
 * storage format, and B is of size [2 x 2] also in column
 * storage format.
 * 
 * THIS ONE HAS THE WORST ACCES PATTERN POSSIBLE. 
 * POSSIBKE FIX: PUT A AND B INTO ROW MAJOR FORMAT FIRST.
 * 
 */
void myMPC_LA_DENSE_MMT2_2_3_2(myMPC_FLOAT *A, myMPC_FLOAT *B, myMPC_FLOAT *L)
{
    int i, j, k, ii, di;
    myMPC_FLOAT ltemp;
    
    ii = 0; di = 0;
    for( i=0; i<2; i++ ){        
        for( j=0; j<=i; j++ ){
            ltemp = 0; 
            for( k=0; k<3; k++ ){
                ltemp += A[k*2+i]*A[k*2+j];
            }			
			for( k=0; k<2; k++ ){
                ltemp += B[k*2+i]*B[k*2+j];
            }
            L[ii+j] = ltemp;
        }
        ii += ++di;
    }
}


/* 
 * Computes r = b - A*x - B*u
 * where A an B are stored in column major format
 */
void myMPC_LA_DENSE_2MVMSUB2_2_3_2(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<2; i++ ){
		r[i] = b[i] - A[k++]*x[0] - B[m++]*u[0];
	}	
	for( j=1; j<3; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] -= A[k++]*x[j];
		}
	}
	
	for( n=1; n<2; n++ ){
		for( i=0; i<2; i++ ){
			r[i] -= B[m++]*u[n];
		}		
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 2.
 */
void myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_FLOAT *L, myMPC_FLOAT *y, myMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    myMPC_FLOAT xel;    
	int start = 1;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 1;
    for( i=1; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 1;
        for( j=1; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }
        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [2 x 2]
 * and stored in column major format. Note the transpose of M!
 */
void myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<2; i++ ){
		r[i] = b[i];
		for( j=0; j<2; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/*
 * Matrix vector multiplication y = b - M'*x where M is of size [2 x 4]
 * and stored in column major format. Note the transpose of M!
 */
void myMPC_LA_DENSE_MTVMSUB_2_4(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *b, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0; 
	for( i=0; i<4; i++ ){
		r[i] = b[i];
		for( j=0; j<2; j++ ){
			r[i] -= A[k++]*x[j];
		}
	}
}


/**
 * Backward Substitution to solve L^T*x = y where L is a
 * lower triangular matrix in triangular storage format.
 * 
 * All involved dimensions are 4.
 */
void myMPC_LA_DENSE_BACKWARDSUB_4(myMPC_FLOAT *L, myMPC_FLOAT *y, myMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    myMPC_FLOAT xel;    
	int start = 6;
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 3;
    for( i=3; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 3;
        for( j=3; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }
        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/*
 * Vector subtraction z = -x - y for vectors of length 29.
 */
void myMPC_LA_VSUB2_29(myMPC_FLOAT *x, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<29; i++){
		z[i] = -x[i] - y[i];
	}
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * lower triangular matrix of size 3 in lower triangular
 * storage format.
 */
void myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_FLOAT *L, myMPC_FLOAT *b, myMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    myMPC_FLOAT y[3];
    myMPC_FLOAT yel,xel;
	int start = 3;
            
    /* first solve Ly = b by forward substitution */
     ii = 0; di = 0;
    for( i=0; i<3; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }
        y[i] = yel / L[ii+i];
        ii += ++di;
    }
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 2;
    for( i=2; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 2;
        for( j=2; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }
        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/**
 * Forward-Backward-Substitution to solve L*L^T*x = b where L is a
 * lower triangular matrix of size 2 in lower triangular
 * storage format.
 */
void myMPC_LA_DENSE_FORWARDBACKWARDSUB_2(myMPC_FLOAT *L, myMPC_FLOAT *b, myMPC_FLOAT *x)
{
    int i, ii, di, j, jj, dj;
    myMPC_FLOAT y[2];
    myMPC_FLOAT yel,xel;
	int start = 1;
            
    /* first solve Ly = b by forward substitution */
     ii = 0; di = 0;
    for( i=0; i<2; i++ ){
        yel = b[i];        
        for( j=0; j<i; j++ ){
            yel -= y[j]*L[ii+j];
        }
        y[i] = yel / L[ii+i];
        ii += ++di;
    }
    
    /* now solve L^T*x = y by backward substitution */
    ii = start; di = 1;
    for( i=1; i>=0; i-- ){        
        xel = y[i];        
        jj = start; dj = 1;
        for( j=1; j>i; j-- ){
            xel -= x[j]*L[jj+i];
            jj -= dj--;
        }
        x[i] = xel / L[ii+i];
        ii -= di--;
    }
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 3,
 * and x has length 3 and is indexed through yidx.
 */
void myMPC_LA_VSUB_INDEXED_3(myMPC_FLOAT *x, int* xidx, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 3.
 */
void myMPC_LA_VSUB3_3(myMPC_FLOAT *u, myMPC_FLOAT *v, myMPC_FLOAT *w, myMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 3
 * and z, x and yidx are of length 3.
 */
void myMPC_LA_VSUB2_INDEXED_3(myMPC_FLOAT *x, myMPC_FLOAT *y, int* yidx, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/*
 * Vector subtraction z = x(xidx) - y where y, z and xidx are of length 2,
 * and x has length 2 and is indexed through yidx.
 */
void myMPC_LA_VSUB_INDEXED_2(myMPC_FLOAT *x, int* xidx, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = x[xidx[i]] - y[i];
	}
}


/*
 * Vector subtraction x = -u.*v - w for vectors of length 2.
 */
void myMPC_LA_VSUB3_2(myMPC_FLOAT *u, myMPC_FLOAT *v, myMPC_FLOAT *w, myMPC_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++){
		x[i] = -u[i]*v[i] - w[i];
	}
}


/*
 * Vector subtraction z = -x - y(yidx) where y is of length 2
 * and z, x and yidx are of length 2.
 */
void myMPC_LA_VSUB2_INDEXED_2(myMPC_FLOAT *x, myMPC_FLOAT *y, int* yidx, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++){
		z[i] = -x[i] - y[yidx[i]];
	}
}


/**
 * Backtracking line search.
 * 
 * First determine the maximum line length by a feasibility line
 * search, i.e. a ~= argmax{ a \in [0...1] s.t. l+a*dl >= 0 and s+a*ds >= 0}.
 *
 * The function returns either the number of iterations or exits the error code
 * myMPC_NOPROGRESS (should be negative).
 */
int myMPC_LINESEARCH_BACKTRACKING_AFFINE(myMPC_FLOAT *l, myMPC_FLOAT *s, myMPC_FLOAT *dl, myMPC_FLOAT *ds, myMPC_FLOAT *a, myMPC_FLOAT *mu_aff)
{
    int i;
	int lsIt=1;    
    myMPC_FLOAT dltemp;
    myMPC_FLOAT dstemp;
    myMPC_FLOAT mya = 1.0;
    myMPC_FLOAT mymu;
        
    while( 1 ){                        

        /* 
         * Compute both snew and wnew together.
         * We compute also mu_affine along the way here, as the
         * values might be in registers, so it should be cheaper.
         */
        mymu = 0;
        for( i=0; i<58; i++ ){
            dltemp = l[i] + mya*dl[i];
            dstemp = s[i] + mya*ds[i];
            if( dltemp < 0 || dstemp < 0 ){
                lsIt++;
                break;
            } else {                
                mymu += dstemp*dltemp;
            }
        }
        
        /* 
         * If no early termination of the for-loop above occurred, we
         * found the required value of a and we can quit the while loop.
         */
        if( i == 58 ){
            break;
        } else {
            mya *= myMPC_SET_LS_SCALE_AFF;
            if( mya < myMPC_SET_LS_MINSTEP ){
                return myMPC_NOPROGRESS;
            }
        }
    }
    
    /* return new values and iteration counter */
    *a = mya;
    *mu_aff = mymu / (myMPC_FLOAT)58;
    return lsIt;
}


/*
 * Vector subtraction x = u.*v - a where a is a scalar
*  and x,u,v are vectors of length 58.
 */
void myMPC_LA_VSUB5_58(myMPC_FLOAT *u, myMPC_FLOAT *v, myMPC_FLOAT a, myMPC_FLOAT *x)
{
	int i;
	for( i=0; i<58; i++){
		x[i] = u[i]*v[i] - a;
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 3,
 * u, su, uidx are of length 3 and v, sv, vidx are of length 3.
 */
void myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_FLOAT *u, myMPC_FLOAT *su, int* uidx, myMPC_FLOAT *v, myMPC_FLOAT *sv, int* vidx, myMPC_FLOAT *x)
{
	int i;
	for( i=0; i<3; i++ ){
		x[i] = 0;
	}
	for( i=0; i<3; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<3; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void myMPC_LA_DENSE_2MVMADD_4_3_3(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<4; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<3; j++ ){		
		for( i=0; i<4; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<3; n++ ){
		for( i=0; i<4; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<2; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<3; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<3; n++ ){
		for( i=0; i<2; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Computes x=0; x(uidx) += u/su; x(vidx) -= v/sv where x is of length 2,
 * u, su, uidx are of length 2 and v, sv, vidx are of length 2.
 */
void myMPC_LA_VSUB6_INDEXED_2_2_2(myMPC_FLOAT *u, myMPC_FLOAT *su, int* uidx, myMPC_FLOAT *v, myMPC_FLOAT *sv, int* vidx, myMPC_FLOAT *x)
{
	int i;
	for( i=0; i<2; i++ ){
		x[i] = 0;
	}
	for( i=0; i<2; i++){
		x[uidx[i]] += u[i]/su[i];
	}
	for( i=0; i<2; i++){
		x[vidx[i]] -= v[i]/sv[i];
	}
}


/* 
 * Computes r = A*x + B*u
 * where A an B are stored in column major format
 */
void myMPC_LA_DENSE_2MVMADD_2_3_2(myMPC_FLOAT *A, myMPC_FLOAT *x, myMPC_FLOAT *B, myMPC_FLOAT *u, myMPC_FLOAT *r)
{
	int i;
	int j;
	int k = 0;
	int m = 0;
	int n;

	for( i=0; i<2; i++ ){
		r[i] = A[k++]*x[0] + B[m++]*u[0];
	}	

	for( j=1; j<3; j++ ){		
		for( i=0; i<2; i++ ){
			r[i] += A[k++]*x[j];
		}
	}
	
	for( n=1; n<2; n++ ){
		for( i=0; i<2; i++ ){
			r[i] += B[m++]*u[n];
		}		
	}
}


/*
 * Vector subtraction z = x - y for vectors of length 29.
 */
void myMPC_LA_VSUB_29(myMPC_FLOAT *x, myMPC_FLOAT *y, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<29; i++){
		z[i] = x[i] - y[i];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_FLOAT *r, myMPC_FLOAT *s, myMPC_FLOAT *u, myMPC_FLOAT *y, int* yidx, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 3 (length of y >= 3).
 */
void myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_FLOAT *r, myMPC_FLOAT *s, myMPC_FLOAT *u, myMPC_FLOAT *y, int* yidx, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<3; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s - u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(myMPC_FLOAT *r, myMPC_FLOAT *s, myMPC_FLOAT *u, myMPC_FLOAT *y, int* yidx, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] - u[i]*y[yidx[i]];
	}
}


/** 
 * Computes z = -r./s + u.*y(y)
 * where all vectors except of y are of length 2 (length of y >= 2).
 */
void myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(myMPC_FLOAT *r, myMPC_FLOAT *s, myMPC_FLOAT *u, myMPC_FLOAT *y, int* yidx, myMPC_FLOAT *z)
{
	int i;
	for( i=0; i<2; i++ ){
		z[i] = -r[i]/s[i] + u[i]*y[yidx[i]];
	}
}


/*
 * Computes ds = -l.\(r + s.*dl) for vectors of length 58.
 */
void myMPC_LA_VSUB7_58(myMPC_FLOAT *l, myMPC_FLOAT *r, myMPC_FLOAT *s, myMPC_FLOAT *dl, myMPC_FLOAT *ds)
{
	int i;
	for( i=0; i<58; i++){
		ds[i] = -(r[i] + s[i]*dl[i])/l[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 29.
 */
void myMPC_LA_VADD_29(myMPC_FLOAT *x, myMPC_FLOAT *y)
{
	int i;
	for( i=0; i<29; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 20.
 */
void myMPC_LA_VADD_20(myMPC_FLOAT *x, myMPC_FLOAT *y)
{
	int i;
	for( i=0; i<20; i++){
		x[i] += y[i];
	}
}


/*
 * Vector addition x = x + y for vectors of length 58.
 */
void myMPC_LA_VADD_58(myMPC_FLOAT *x, myMPC_FLOAT *y)
{
	int i;
	for( i=0; i<58; i++){
		x[i] += y[i];
	}
}


/**
 * Backtracking line search for combined predictor/corrector step.
 * Update on variables with safety factor gamma (to keep us away from
 * boundary).
 */
int myMPC_LINESEARCH_BACKTRACKING_COMBINED(myMPC_FLOAT *z, myMPC_FLOAT *v, myMPC_FLOAT *l, myMPC_FLOAT *s, myMPC_FLOAT *dz, myMPC_FLOAT *dv, myMPC_FLOAT *dl, myMPC_FLOAT *ds, myMPC_FLOAT *a, myMPC_FLOAT *mu)
{
    int i, lsIt=1;       
    myMPC_FLOAT dltemp;
    myMPC_FLOAT dstemp;    
    myMPC_FLOAT a_gamma;
            
    *a = 1.0;
    while( 1 ){                        

        /* check whether search criterion is fulfilled */
        for( i=0; i<58; i++ ){
            dltemp = l[i] + (*a)*dl[i];
            dstemp = s[i] + (*a)*ds[i];
            if( dltemp < 0 || dstemp < 0 ){
                lsIt++;
                break;
            }
        }
        
        /* 
         * If no early termination of the for-loop above occurred, we
         * found the required value of a and we can quit the while loop.
         */
        if( i == 58 ){
            break;
        } else {
            *a *= myMPC_SET_LS_SCALE;
            if( *a < myMPC_SET_LS_MINSTEP ){
                return myMPC_NOPROGRESS;
            }
        }
    }
    
    /* update variables with safety margin */
    a_gamma = (*a)*myMPC_SET_LS_MAXSTEP;
    
    /* primal variables */
    for( i=0; i<29; i++ ){
        z[i] += a_gamma*dz[i];
    }
    
    /* equality constraint multipliers */
    for( i=0; i<20; i++ ){
        v[i] += a_gamma*dv[i];
    }
    
    /* inequality constraint multipliers & slacks, also update mu */
    *mu = 0;
    for( i=0; i<58; i++ ){
        dltemp = l[i] + a_gamma*dl[i]; l[i] = dltemp;
        dstemp = s[i] + a_gamma*ds[i]; s[i] = dstemp;
        *mu += dltemp*dstemp;
    }
    
    *a = a_gamma;
    *mu /= (myMPC_FLOAT)58;
    return lsIt;
}




/* VARIABLE DEFINITIONS ------------------------------------------------ */
myMPC_FLOAT myMPC_z[29];
myMPC_FLOAT myMPC_v[20];
myMPC_FLOAT myMPC_dz_aff[29];
myMPC_FLOAT myMPC_dv_aff[20];
myMPC_FLOAT myMPC_grad_cost[29];
myMPC_FLOAT myMPC_grad_eq[29];
myMPC_FLOAT myMPC_rd[29];
myMPC_FLOAT myMPC_l[58];
myMPC_FLOAT myMPC_s[58];
myMPC_FLOAT myMPC_lbys[58];
myMPC_FLOAT myMPC_dl_aff[58];
myMPC_FLOAT myMPC_ds_aff[58];
myMPC_FLOAT myMPC_dz_cc[29];
myMPC_FLOAT myMPC_dv_cc[20];
myMPC_FLOAT myMPC_dl_cc[58];
myMPC_FLOAT myMPC_ds_cc[58];
myMPC_FLOAT myMPC_ccrhs[58];
myMPC_FLOAT myMPC_grad_ineq[29];
myMPC_FLOAT myMPC_H0[9] = {1.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000};
myMPC_FLOAT myMPC_f0[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z0 = myMPC_z + 0;
myMPC_FLOAT* myMPC_dzaff0 = myMPC_dz_aff + 0;
myMPC_FLOAT* myMPC_dzcc0 = myMPC_dz_cc + 0;
myMPC_FLOAT* myMPC_rd0 = myMPC_rd + 0;
myMPC_FLOAT myMPC_Lbyrd0[3];
myMPC_FLOAT* myMPC_grad_cost0 = myMPC_grad_cost + 0;
myMPC_FLOAT* myMPC_grad_eq0 = myMPC_grad_eq + 0;
myMPC_FLOAT* myMPC_grad_ineq0 = myMPC_grad_ineq + 0;
myMPC_FLOAT myMPC_ctv0[3];
myMPC_FLOAT myMPC_Phi0[6];
myMPC_FLOAT myMPC_C0[12] = {1.0000000000000000E+000, 0.0000000000000000E+000, 1.1000000000000001E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 1.0000000000000000E+000, 1.0000000000000000E+000, 1.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 1.0000000000000000E+000, 5.0000000000000000E-001};
myMPC_FLOAT* myMPC_v0 = myMPC_v + 0;
myMPC_FLOAT myMPC_re0[4];
myMPC_FLOAT myMPC_beta0[4];
myMPC_FLOAT myMPC_betacc0[4];
myMPC_FLOAT* myMPC_dvaff0 = myMPC_dv_aff + 0;
myMPC_FLOAT* myMPC_dvcc0 = myMPC_dv_cc + 0;
myMPC_FLOAT myMPC_V0[12];
myMPC_FLOAT myMPC_Yd0[10];
myMPC_FLOAT myMPC_Ld0[10];
myMPC_FLOAT myMPC_yy0[4];
myMPC_FLOAT myMPC_bmy0[4];
myMPC_FLOAT myMPC_lb0[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx0[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb0 = myMPC_l + 0;
myMPC_FLOAT* myMPC_slb0 = myMPC_s + 0;
myMPC_FLOAT* myMPC_llbbyslb0 = myMPC_lbys + 0;
myMPC_FLOAT myMPC_rilb0[3];
myMPC_FLOAT* myMPC_dllbaff0 = myMPC_dl_aff + 0;
myMPC_FLOAT* myMPC_dslbaff0 = myMPC_ds_aff + 0;
myMPC_FLOAT* myMPC_dllbcc0 = myMPC_dl_cc + 0;
myMPC_FLOAT* myMPC_dslbcc0 = myMPC_ds_cc + 0;
myMPC_FLOAT* myMPC_ccrhsl0 = myMPC_ccrhs + 0;
myMPC_FLOAT myMPC_ub0[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx0[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub0 = myMPC_l + 3;
myMPC_FLOAT* myMPC_sub0 = myMPC_s + 3;
myMPC_FLOAT* myMPC_lubbysub0 = myMPC_lbys + 3;
myMPC_FLOAT myMPC_riub0[3];
myMPC_FLOAT* myMPC_dlubaff0 = myMPC_dl_aff + 3;
myMPC_FLOAT* myMPC_dsubaff0 = myMPC_ds_aff + 3;
myMPC_FLOAT* myMPC_dlubcc0 = myMPC_dl_cc + 3;
myMPC_FLOAT* myMPC_dsubcc0 = myMPC_ds_cc + 3;
myMPC_FLOAT* myMPC_ccrhsub0 = myMPC_ccrhs + 3;
myMPC_FLOAT myMPC_f1[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z1 = myMPC_z + 3;
myMPC_FLOAT* myMPC_dzaff1 = myMPC_dz_aff + 3;
myMPC_FLOAT* myMPC_dzcc1 = myMPC_dz_cc + 3;
myMPC_FLOAT* myMPC_rd1 = myMPC_rd + 3;
myMPC_FLOAT myMPC_Lbyrd1[3];
myMPC_FLOAT* myMPC_grad_cost1 = myMPC_grad_cost + 3;
myMPC_FLOAT* myMPC_grad_eq1 = myMPC_grad_eq + 3;
myMPC_FLOAT* myMPC_grad_ineq1 = myMPC_grad_ineq + 3;
myMPC_FLOAT myMPC_ctv1[3];
myMPC_FLOAT myMPC_Phi1[6];
myMPC_FLOAT myMPC_C1[6] = {1.1000000000000001E+000, 0.0000000000000000E+000, 
1.0000000000000000E+000, 1.0000000000000000E+000, 
1.0000000000000000E+000, 5.0000000000000000E-001};
myMPC_FLOAT* myMPC_v1 = myMPC_v + 4;
myMPC_FLOAT myMPC_re1[2];
myMPC_FLOAT myMPC_beta1[2];
myMPC_FLOAT myMPC_betacc1[2];
myMPC_FLOAT* myMPC_dvaff1 = myMPC_dv_aff + 4;
myMPC_FLOAT* myMPC_dvcc1 = myMPC_dv_cc + 4;
myMPC_FLOAT myMPC_V1[6];
myMPC_FLOAT myMPC_Yd1[3];
myMPC_FLOAT myMPC_Ld1[3];
myMPC_FLOAT myMPC_yy1[2];
myMPC_FLOAT myMPC_bmy1[2];
myMPC_FLOAT myMPC_c1[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_D1[12] = {0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, -1.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W1[12];
myMPC_FLOAT myMPC_lb1[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx1[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb1 = myMPC_l + 6;
myMPC_FLOAT* myMPC_slb1 = myMPC_s + 6;
myMPC_FLOAT* myMPC_llbbyslb1 = myMPC_lbys + 6;
myMPC_FLOAT myMPC_rilb1[3];
myMPC_FLOAT* myMPC_dllbaff1 = myMPC_dl_aff + 6;
myMPC_FLOAT* myMPC_dslbaff1 = myMPC_ds_aff + 6;
myMPC_FLOAT* myMPC_dllbcc1 = myMPC_dl_cc + 6;
myMPC_FLOAT* myMPC_dslbcc1 = myMPC_ds_cc + 6;
myMPC_FLOAT* myMPC_ccrhsl1 = myMPC_ccrhs + 6;
myMPC_FLOAT myMPC_ub1[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx1[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub1 = myMPC_l + 9;
myMPC_FLOAT* myMPC_sub1 = myMPC_s + 9;
myMPC_FLOAT* myMPC_lubbysub1 = myMPC_lbys + 9;
myMPC_FLOAT myMPC_riub1[3];
myMPC_FLOAT* myMPC_dlubaff1 = myMPC_dl_aff + 9;
myMPC_FLOAT* myMPC_dsubaff1 = myMPC_ds_aff + 9;
myMPC_FLOAT* myMPC_dlubcc1 = myMPC_dl_cc + 9;
myMPC_FLOAT* myMPC_dsubcc1 = myMPC_ds_cc + 9;
myMPC_FLOAT* myMPC_ccrhsub1 = myMPC_ccrhs + 9;
myMPC_FLOAT myMPC_Ysd1[8];
myMPC_FLOAT myMPC_Lsd1[8];
myMPC_FLOAT myMPC_f2[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z2 = myMPC_z + 6;
myMPC_FLOAT* myMPC_dzaff2 = myMPC_dz_aff + 6;
myMPC_FLOAT* myMPC_dzcc2 = myMPC_dz_cc + 6;
myMPC_FLOAT* myMPC_rd2 = myMPC_rd + 6;
myMPC_FLOAT myMPC_Lbyrd2[3];
myMPC_FLOAT* myMPC_grad_cost2 = myMPC_grad_cost + 6;
myMPC_FLOAT* myMPC_grad_eq2 = myMPC_grad_eq + 6;
myMPC_FLOAT* myMPC_grad_ineq2 = myMPC_grad_ineq + 6;
myMPC_FLOAT myMPC_ctv2[3];
myMPC_FLOAT myMPC_Phi2[6];
myMPC_FLOAT* myMPC_v2 = myMPC_v + 6;
myMPC_FLOAT myMPC_re2[2];
myMPC_FLOAT myMPC_beta2[2];
myMPC_FLOAT myMPC_betacc2[2];
myMPC_FLOAT* myMPC_dvaff2 = myMPC_dv_aff + 6;
myMPC_FLOAT* myMPC_dvcc2 = myMPC_dv_cc + 6;
myMPC_FLOAT myMPC_V2[6];
myMPC_FLOAT myMPC_Yd2[3];
myMPC_FLOAT myMPC_Ld2[3];
myMPC_FLOAT myMPC_yy2[2];
myMPC_FLOAT myMPC_bmy2[2];
myMPC_FLOAT myMPC_c2[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_D2[6] = {-1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, -1.0000000000000000E+000, 
0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W2[6];
myMPC_FLOAT myMPC_lb2[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx2[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb2 = myMPC_l + 12;
myMPC_FLOAT* myMPC_slb2 = myMPC_s + 12;
myMPC_FLOAT* myMPC_llbbyslb2 = myMPC_lbys + 12;
myMPC_FLOAT myMPC_rilb2[3];
myMPC_FLOAT* myMPC_dllbaff2 = myMPC_dl_aff + 12;
myMPC_FLOAT* myMPC_dslbaff2 = myMPC_ds_aff + 12;
myMPC_FLOAT* myMPC_dllbcc2 = myMPC_dl_cc + 12;
myMPC_FLOAT* myMPC_dslbcc2 = myMPC_ds_cc + 12;
myMPC_FLOAT* myMPC_ccrhsl2 = myMPC_ccrhs + 12;
myMPC_FLOAT myMPC_ub2[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx2[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub2 = myMPC_l + 15;
myMPC_FLOAT* myMPC_sub2 = myMPC_s + 15;
myMPC_FLOAT* myMPC_lubbysub2 = myMPC_lbys + 15;
myMPC_FLOAT myMPC_riub2[3];
myMPC_FLOAT* myMPC_dlubaff2 = myMPC_dl_aff + 15;
myMPC_FLOAT* myMPC_dsubaff2 = myMPC_ds_aff + 15;
myMPC_FLOAT* myMPC_dlubcc2 = myMPC_dl_cc + 15;
myMPC_FLOAT* myMPC_dsubcc2 = myMPC_ds_cc + 15;
myMPC_FLOAT* myMPC_ccrhsub2 = myMPC_ccrhs + 15;
myMPC_FLOAT myMPC_Ysd2[4];
myMPC_FLOAT myMPC_Lsd2[4];
myMPC_FLOAT myMPC_f3[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z3 = myMPC_z + 9;
myMPC_FLOAT* myMPC_dzaff3 = myMPC_dz_aff + 9;
myMPC_FLOAT* myMPC_dzcc3 = myMPC_dz_cc + 9;
myMPC_FLOAT* myMPC_rd3 = myMPC_rd + 9;
myMPC_FLOAT myMPC_Lbyrd3[3];
myMPC_FLOAT* myMPC_grad_cost3 = myMPC_grad_cost + 9;
myMPC_FLOAT* myMPC_grad_eq3 = myMPC_grad_eq + 9;
myMPC_FLOAT* myMPC_grad_ineq3 = myMPC_grad_ineq + 9;
myMPC_FLOAT myMPC_ctv3[3];
myMPC_FLOAT myMPC_Phi3[6];
myMPC_FLOAT* myMPC_v3 = myMPC_v + 8;
myMPC_FLOAT myMPC_re3[2];
myMPC_FLOAT myMPC_beta3[2];
myMPC_FLOAT myMPC_betacc3[2];
myMPC_FLOAT* myMPC_dvaff3 = myMPC_dv_aff + 8;
myMPC_FLOAT* myMPC_dvcc3 = myMPC_dv_cc + 8;
myMPC_FLOAT myMPC_V3[6];
myMPC_FLOAT myMPC_Yd3[3];
myMPC_FLOAT myMPC_Ld3[3];
myMPC_FLOAT myMPC_yy3[2];
myMPC_FLOAT myMPC_bmy3[2];
myMPC_FLOAT myMPC_c3[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W3[6];
myMPC_FLOAT myMPC_lb3[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx3[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb3 = myMPC_l + 18;
myMPC_FLOAT* myMPC_slb3 = myMPC_s + 18;
myMPC_FLOAT* myMPC_llbbyslb3 = myMPC_lbys + 18;
myMPC_FLOAT myMPC_rilb3[3];
myMPC_FLOAT* myMPC_dllbaff3 = myMPC_dl_aff + 18;
myMPC_FLOAT* myMPC_dslbaff3 = myMPC_ds_aff + 18;
myMPC_FLOAT* myMPC_dllbcc3 = myMPC_dl_cc + 18;
myMPC_FLOAT* myMPC_dslbcc3 = myMPC_ds_cc + 18;
myMPC_FLOAT* myMPC_ccrhsl3 = myMPC_ccrhs + 18;
myMPC_FLOAT myMPC_ub3[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx3[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub3 = myMPC_l + 21;
myMPC_FLOAT* myMPC_sub3 = myMPC_s + 21;
myMPC_FLOAT* myMPC_lubbysub3 = myMPC_lbys + 21;
myMPC_FLOAT myMPC_riub3[3];
myMPC_FLOAT* myMPC_dlubaff3 = myMPC_dl_aff + 21;
myMPC_FLOAT* myMPC_dsubaff3 = myMPC_ds_aff + 21;
myMPC_FLOAT* myMPC_dlubcc3 = myMPC_dl_cc + 21;
myMPC_FLOAT* myMPC_dsubcc3 = myMPC_ds_cc + 21;
myMPC_FLOAT* myMPC_ccrhsub3 = myMPC_ccrhs + 21;
myMPC_FLOAT myMPC_Ysd3[4];
myMPC_FLOAT myMPC_Lsd3[4];
myMPC_FLOAT myMPC_f4[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z4 = myMPC_z + 12;
myMPC_FLOAT* myMPC_dzaff4 = myMPC_dz_aff + 12;
myMPC_FLOAT* myMPC_dzcc4 = myMPC_dz_cc + 12;
myMPC_FLOAT* myMPC_rd4 = myMPC_rd + 12;
myMPC_FLOAT myMPC_Lbyrd4[3];
myMPC_FLOAT* myMPC_grad_cost4 = myMPC_grad_cost + 12;
myMPC_FLOAT* myMPC_grad_eq4 = myMPC_grad_eq + 12;
myMPC_FLOAT* myMPC_grad_ineq4 = myMPC_grad_ineq + 12;
myMPC_FLOAT myMPC_ctv4[3];
myMPC_FLOAT myMPC_Phi4[6];
myMPC_FLOAT* myMPC_v4 = myMPC_v + 10;
myMPC_FLOAT myMPC_re4[2];
myMPC_FLOAT myMPC_beta4[2];
myMPC_FLOAT myMPC_betacc4[2];
myMPC_FLOAT* myMPC_dvaff4 = myMPC_dv_aff + 10;
myMPC_FLOAT* myMPC_dvcc4 = myMPC_dv_cc + 10;
myMPC_FLOAT myMPC_V4[6];
myMPC_FLOAT myMPC_Yd4[3];
myMPC_FLOAT myMPC_Ld4[3];
myMPC_FLOAT myMPC_yy4[2];
myMPC_FLOAT myMPC_bmy4[2];
myMPC_FLOAT myMPC_c4[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W4[6];
myMPC_FLOAT myMPC_lb4[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx4[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb4 = myMPC_l + 24;
myMPC_FLOAT* myMPC_slb4 = myMPC_s + 24;
myMPC_FLOAT* myMPC_llbbyslb4 = myMPC_lbys + 24;
myMPC_FLOAT myMPC_rilb4[3];
myMPC_FLOAT* myMPC_dllbaff4 = myMPC_dl_aff + 24;
myMPC_FLOAT* myMPC_dslbaff4 = myMPC_ds_aff + 24;
myMPC_FLOAT* myMPC_dllbcc4 = myMPC_dl_cc + 24;
myMPC_FLOAT* myMPC_dslbcc4 = myMPC_ds_cc + 24;
myMPC_FLOAT* myMPC_ccrhsl4 = myMPC_ccrhs + 24;
myMPC_FLOAT myMPC_ub4[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx4[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub4 = myMPC_l + 27;
myMPC_FLOAT* myMPC_sub4 = myMPC_s + 27;
myMPC_FLOAT* myMPC_lubbysub4 = myMPC_lbys + 27;
myMPC_FLOAT myMPC_riub4[3];
myMPC_FLOAT* myMPC_dlubaff4 = myMPC_dl_aff + 27;
myMPC_FLOAT* myMPC_dsubaff4 = myMPC_ds_aff + 27;
myMPC_FLOAT* myMPC_dlubcc4 = myMPC_dl_cc + 27;
myMPC_FLOAT* myMPC_dsubcc4 = myMPC_ds_cc + 27;
myMPC_FLOAT* myMPC_ccrhsub4 = myMPC_ccrhs + 27;
myMPC_FLOAT myMPC_Ysd4[4];
myMPC_FLOAT myMPC_Lsd4[4];
myMPC_FLOAT myMPC_f5[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z5 = myMPC_z + 15;
myMPC_FLOAT* myMPC_dzaff5 = myMPC_dz_aff + 15;
myMPC_FLOAT* myMPC_dzcc5 = myMPC_dz_cc + 15;
myMPC_FLOAT* myMPC_rd5 = myMPC_rd + 15;
myMPC_FLOAT myMPC_Lbyrd5[3];
myMPC_FLOAT* myMPC_grad_cost5 = myMPC_grad_cost + 15;
myMPC_FLOAT* myMPC_grad_eq5 = myMPC_grad_eq + 15;
myMPC_FLOAT* myMPC_grad_ineq5 = myMPC_grad_ineq + 15;
myMPC_FLOAT myMPC_ctv5[3];
myMPC_FLOAT myMPC_Phi5[6];
myMPC_FLOAT* myMPC_v5 = myMPC_v + 12;
myMPC_FLOAT myMPC_re5[2];
myMPC_FLOAT myMPC_beta5[2];
myMPC_FLOAT myMPC_betacc5[2];
myMPC_FLOAT* myMPC_dvaff5 = myMPC_dv_aff + 12;
myMPC_FLOAT* myMPC_dvcc5 = myMPC_dv_cc + 12;
myMPC_FLOAT myMPC_V5[6];
myMPC_FLOAT myMPC_Yd5[3];
myMPC_FLOAT myMPC_Ld5[3];
myMPC_FLOAT myMPC_yy5[2];
myMPC_FLOAT myMPC_bmy5[2];
myMPC_FLOAT myMPC_c5[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W5[6];
myMPC_FLOAT myMPC_lb5[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx5[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb5 = myMPC_l + 30;
myMPC_FLOAT* myMPC_slb5 = myMPC_s + 30;
myMPC_FLOAT* myMPC_llbbyslb5 = myMPC_lbys + 30;
myMPC_FLOAT myMPC_rilb5[3];
myMPC_FLOAT* myMPC_dllbaff5 = myMPC_dl_aff + 30;
myMPC_FLOAT* myMPC_dslbaff5 = myMPC_ds_aff + 30;
myMPC_FLOAT* myMPC_dllbcc5 = myMPC_dl_cc + 30;
myMPC_FLOAT* myMPC_dslbcc5 = myMPC_ds_cc + 30;
myMPC_FLOAT* myMPC_ccrhsl5 = myMPC_ccrhs + 30;
myMPC_FLOAT myMPC_ub5[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx5[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub5 = myMPC_l + 33;
myMPC_FLOAT* myMPC_sub5 = myMPC_s + 33;
myMPC_FLOAT* myMPC_lubbysub5 = myMPC_lbys + 33;
myMPC_FLOAT myMPC_riub5[3];
myMPC_FLOAT* myMPC_dlubaff5 = myMPC_dl_aff + 33;
myMPC_FLOAT* myMPC_dsubaff5 = myMPC_ds_aff + 33;
myMPC_FLOAT* myMPC_dlubcc5 = myMPC_dl_cc + 33;
myMPC_FLOAT* myMPC_dsubcc5 = myMPC_ds_cc + 33;
myMPC_FLOAT* myMPC_ccrhsub5 = myMPC_ccrhs + 33;
myMPC_FLOAT myMPC_Ysd5[4];
myMPC_FLOAT myMPC_Lsd5[4];
myMPC_FLOAT myMPC_f6[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z6 = myMPC_z + 18;
myMPC_FLOAT* myMPC_dzaff6 = myMPC_dz_aff + 18;
myMPC_FLOAT* myMPC_dzcc6 = myMPC_dz_cc + 18;
myMPC_FLOAT* myMPC_rd6 = myMPC_rd + 18;
myMPC_FLOAT myMPC_Lbyrd6[3];
myMPC_FLOAT* myMPC_grad_cost6 = myMPC_grad_cost + 18;
myMPC_FLOAT* myMPC_grad_eq6 = myMPC_grad_eq + 18;
myMPC_FLOAT* myMPC_grad_ineq6 = myMPC_grad_ineq + 18;
myMPC_FLOAT myMPC_ctv6[3];
myMPC_FLOAT myMPC_Phi6[6];
myMPC_FLOAT* myMPC_v6 = myMPC_v + 14;
myMPC_FLOAT myMPC_re6[2];
myMPC_FLOAT myMPC_beta6[2];
myMPC_FLOAT myMPC_betacc6[2];
myMPC_FLOAT* myMPC_dvaff6 = myMPC_dv_aff + 14;
myMPC_FLOAT* myMPC_dvcc6 = myMPC_dv_cc + 14;
myMPC_FLOAT myMPC_V6[6];
myMPC_FLOAT myMPC_Yd6[3];
myMPC_FLOAT myMPC_Ld6[3];
myMPC_FLOAT myMPC_yy6[2];
myMPC_FLOAT myMPC_bmy6[2];
myMPC_FLOAT myMPC_c6[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W6[6];
myMPC_FLOAT myMPC_lb6[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx6[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb6 = myMPC_l + 36;
myMPC_FLOAT* myMPC_slb6 = myMPC_s + 36;
myMPC_FLOAT* myMPC_llbbyslb6 = myMPC_lbys + 36;
myMPC_FLOAT myMPC_rilb6[3];
myMPC_FLOAT* myMPC_dllbaff6 = myMPC_dl_aff + 36;
myMPC_FLOAT* myMPC_dslbaff6 = myMPC_ds_aff + 36;
myMPC_FLOAT* myMPC_dllbcc6 = myMPC_dl_cc + 36;
myMPC_FLOAT* myMPC_dslbcc6 = myMPC_ds_cc + 36;
myMPC_FLOAT* myMPC_ccrhsl6 = myMPC_ccrhs + 36;
myMPC_FLOAT myMPC_ub6[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx6[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub6 = myMPC_l + 39;
myMPC_FLOAT* myMPC_sub6 = myMPC_s + 39;
myMPC_FLOAT* myMPC_lubbysub6 = myMPC_lbys + 39;
myMPC_FLOAT myMPC_riub6[3];
myMPC_FLOAT* myMPC_dlubaff6 = myMPC_dl_aff + 39;
myMPC_FLOAT* myMPC_dsubaff6 = myMPC_ds_aff + 39;
myMPC_FLOAT* myMPC_dlubcc6 = myMPC_dl_cc + 39;
myMPC_FLOAT* myMPC_dsubcc6 = myMPC_ds_cc + 39;
myMPC_FLOAT* myMPC_ccrhsub6 = myMPC_ccrhs + 39;
myMPC_FLOAT myMPC_Ysd6[4];
myMPC_FLOAT myMPC_Lsd6[4];
myMPC_FLOAT myMPC_f7[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z7 = myMPC_z + 21;
myMPC_FLOAT* myMPC_dzaff7 = myMPC_dz_aff + 21;
myMPC_FLOAT* myMPC_dzcc7 = myMPC_dz_cc + 21;
myMPC_FLOAT* myMPC_rd7 = myMPC_rd + 21;
myMPC_FLOAT myMPC_Lbyrd7[3];
myMPC_FLOAT* myMPC_grad_cost7 = myMPC_grad_cost + 21;
myMPC_FLOAT* myMPC_grad_eq7 = myMPC_grad_eq + 21;
myMPC_FLOAT* myMPC_grad_ineq7 = myMPC_grad_ineq + 21;
myMPC_FLOAT myMPC_ctv7[3];
myMPC_FLOAT myMPC_Phi7[6];
myMPC_FLOAT* myMPC_v7 = myMPC_v + 16;
myMPC_FLOAT myMPC_re7[2];
myMPC_FLOAT myMPC_beta7[2];
myMPC_FLOAT myMPC_betacc7[2];
myMPC_FLOAT* myMPC_dvaff7 = myMPC_dv_aff + 16;
myMPC_FLOAT* myMPC_dvcc7 = myMPC_dv_cc + 16;
myMPC_FLOAT myMPC_V7[6];
myMPC_FLOAT myMPC_Yd7[3];
myMPC_FLOAT myMPC_Ld7[3];
myMPC_FLOAT myMPC_yy7[2];
myMPC_FLOAT myMPC_bmy7[2];
myMPC_FLOAT myMPC_c7[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W7[6];
myMPC_FLOAT myMPC_lb7[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx7[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb7 = myMPC_l + 42;
myMPC_FLOAT* myMPC_slb7 = myMPC_s + 42;
myMPC_FLOAT* myMPC_llbbyslb7 = myMPC_lbys + 42;
myMPC_FLOAT myMPC_rilb7[3];
myMPC_FLOAT* myMPC_dllbaff7 = myMPC_dl_aff + 42;
myMPC_FLOAT* myMPC_dslbaff7 = myMPC_ds_aff + 42;
myMPC_FLOAT* myMPC_dllbcc7 = myMPC_dl_cc + 42;
myMPC_FLOAT* myMPC_dslbcc7 = myMPC_ds_cc + 42;
myMPC_FLOAT* myMPC_ccrhsl7 = myMPC_ccrhs + 42;
myMPC_FLOAT myMPC_ub7[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx7[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub7 = myMPC_l + 45;
myMPC_FLOAT* myMPC_sub7 = myMPC_s + 45;
myMPC_FLOAT* myMPC_lubbysub7 = myMPC_lbys + 45;
myMPC_FLOAT myMPC_riub7[3];
myMPC_FLOAT* myMPC_dlubaff7 = myMPC_dl_aff + 45;
myMPC_FLOAT* myMPC_dsubaff7 = myMPC_ds_aff + 45;
myMPC_FLOAT* myMPC_dlubcc7 = myMPC_dl_cc + 45;
myMPC_FLOAT* myMPC_dsubcc7 = myMPC_ds_cc + 45;
myMPC_FLOAT* myMPC_ccrhsub7 = myMPC_ccrhs + 45;
myMPC_FLOAT myMPC_Ysd7[4];
myMPC_FLOAT myMPC_Lsd7[4];
myMPC_FLOAT myMPC_f8[3] = {0.0000000000000000E+000, 0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z8 = myMPC_z + 24;
myMPC_FLOAT* myMPC_dzaff8 = myMPC_dz_aff + 24;
myMPC_FLOAT* myMPC_dzcc8 = myMPC_dz_cc + 24;
myMPC_FLOAT* myMPC_rd8 = myMPC_rd + 24;
myMPC_FLOAT myMPC_Lbyrd8[3];
myMPC_FLOAT* myMPC_grad_cost8 = myMPC_grad_cost + 24;
myMPC_FLOAT* myMPC_grad_eq8 = myMPC_grad_eq + 24;
myMPC_FLOAT* myMPC_grad_ineq8 = myMPC_grad_ineq + 24;
myMPC_FLOAT myMPC_ctv8[3];
myMPC_FLOAT myMPC_Phi8[6];
myMPC_FLOAT* myMPC_v8 = myMPC_v + 18;
myMPC_FLOAT myMPC_re8[2];
myMPC_FLOAT myMPC_beta8[2];
myMPC_FLOAT myMPC_betacc8[2];
myMPC_FLOAT* myMPC_dvaff8 = myMPC_dv_aff + 18;
myMPC_FLOAT* myMPC_dvcc8 = myMPC_dv_cc + 18;
myMPC_FLOAT myMPC_V8[6];
myMPC_FLOAT myMPC_Yd8[3];
myMPC_FLOAT myMPC_Ld8[3];
myMPC_FLOAT myMPC_yy8[2];
myMPC_FLOAT myMPC_bmy8[2];
myMPC_FLOAT myMPC_c8[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT myMPC_W8[6];
myMPC_FLOAT myMPC_lb8[3] = {-5.0000000000000000E+000, -5.0000000000000000E+000, -5.0000000000000000E-001};
int myMPC_lbIdx8[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_llb8 = myMPC_l + 48;
myMPC_FLOAT* myMPC_slb8 = myMPC_s + 48;
myMPC_FLOAT* myMPC_llbbyslb8 = myMPC_lbys + 48;
myMPC_FLOAT myMPC_rilb8[3];
myMPC_FLOAT* myMPC_dllbaff8 = myMPC_dl_aff + 48;
myMPC_FLOAT* myMPC_dslbaff8 = myMPC_ds_aff + 48;
myMPC_FLOAT* myMPC_dllbcc8 = myMPC_dl_cc + 48;
myMPC_FLOAT* myMPC_dslbcc8 = myMPC_ds_cc + 48;
myMPC_FLOAT* myMPC_ccrhsl8 = myMPC_ccrhs + 48;
myMPC_FLOAT myMPC_ub8[3] = {5.0000000000000000E+000, 5.0000000000000000E+000, 5.0000000000000000E-001};
int myMPC_ubIdx8[3] = {0, 1, 2};
myMPC_FLOAT* myMPC_lub8 = myMPC_l + 51;
myMPC_FLOAT* myMPC_sub8 = myMPC_s + 51;
myMPC_FLOAT* myMPC_lubbysub8 = myMPC_lbys + 51;
myMPC_FLOAT myMPC_riub8[3];
myMPC_FLOAT* myMPC_dlubaff8 = myMPC_dl_aff + 51;
myMPC_FLOAT* myMPC_dsubaff8 = myMPC_ds_aff + 51;
myMPC_FLOAT* myMPC_dlubcc8 = myMPC_dl_cc + 51;
myMPC_FLOAT* myMPC_dsubcc8 = myMPC_ds_cc + 51;
myMPC_FLOAT* myMPC_ccrhsub8 = myMPC_ccrhs + 51;
myMPC_FLOAT myMPC_Ysd8[4];
myMPC_FLOAT myMPC_Lsd8[4];
myMPC_FLOAT myMPC_H9[4] = {2.0239004884654217E+000, 2.6945484784212870E-001, 
2.6945484784212870E-001, 2.6529099408645656E+000};
myMPC_FLOAT myMPC_f9[2] = {0.0000000000000000E+000, 0.0000000000000000E+000};
myMPC_FLOAT* myMPC_z9 = myMPC_z + 27;
myMPC_FLOAT* myMPC_dzaff9 = myMPC_dz_aff + 27;
myMPC_FLOAT* myMPC_dzcc9 = myMPC_dz_cc + 27;
myMPC_FLOAT* myMPC_rd9 = myMPC_rd + 27;
myMPC_FLOAT myMPC_Lbyrd9[2];
myMPC_FLOAT* myMPC_grad_cost9 = myMPC_grad_cost + 27;
myMPC_FLOAT* myMPC_grad_eq9 = myMPC_grad_eq + 27;
myMPC_FLOAT* myMPC_grad_ineq9 = myMPC_grad_ineq + 27;
myMPC_FLOAT myMPC_ctv9[2];
myMPC_FLOAT myMPC_Phi9[3];
myMPC_FLOAT myMPC_D9[4] = {-1.0000000000000000E+000, 0.0000000000000000E+000, 
0.0000000000000000E+000, -1.0000000000000000E+000};
myMPC_FLOAT myMPC_W9[4];
myMPC_FLOAT myMPC_lb9[2] = {-5.0000000000000000E+000, -5.0000000000000000E+000};
int myMPC_lbIdx9[2] = {0, 1};
myMPC_FLOAT* myMPC_llb9 = myMPC_l + 54;
myMPC_FLOAT* myMPC_slb9 = myMPC_s + 54;
myMPC_FLOAT* myMPC_llbbyslb9 = myMPC_lbys + 54;
myMPC_FLOAT myMPC_rilb9[2];
myMPC_FLOAT* myMPC_dllbaff9 = myMPC_dl_aff + 54;
myMPC_FLOAT* myMPC_dslbaff9 = myMPC_ds_aff + 54;
myMPC_FLOAT* myMPC_dllbcc9 = myMPC_dl_cc + 54;
myMPC_FLOAT* myMPC_dslbcc9 = myMPC_ds_cc + 54;
myMPC_FLOAT* myMPC_ccrhsl9 = myMPC_ccrhs + 54;
myMPC_FLOAT myMPC_ub9[2] = {5.0000000000000000E+000, 5.0000000000000000E+000};
int myMPC_ubIdx9[2] = {0, 1};
myMPC_FLOAT* myMPC_lub9 = myMPC_l + 56;
myMPC_FLOAT* myMPC_sub9 = myMPC_s + 56;
myMPC_FLOAT* myMPC_lubbysub9 = myMPC_lbys + 56;
myMPC_FLOAT myMPC_riub9[2];
myMPC_FLOAT* myMPC_dlubaff9 = myMPC_dl_aff + 56;
myMPC_FLOAT* myMPC_dsubaff9 = myMPC_ds_aff + 56;
myMPC_FLOAT* myMPC_dlubcc9 = myMPC_dl_cc + 56;
myMPC_FLOAT* myMPC_dsubcc9 = myMPC_ds_cc + 56;
myMPC_FLOAT* myMPC_ccrhsub9 = myMPC_ccrhs + 56;
myMPC_FLOAT musigma;
myMPC_FLOAT sigma_3rdroot;


/* SOLVER CODE --------------------------------------------------------- */
int myMPC_solve(myMPC_params* params, myMPC_output* output, myMPC_info* info)
{	
int exitcode;
int i;
#if myMPC_SET_PRINTLEVEL > 0
	myMPC_timer solvertimer;
	myMPC_tic(&solvertimer);
#endif
/* FUNCTION CALLS INTO LA LIBRARY -------------------------------------- */
info->it = 0;
myMPC_LA_INITIALIZEVECTOR_29(myMPC_z, 0);
myMPC_LA_INITIALIZEVECTOR_20(myMPC_v, 1);
myMPC_LA_INITIALIZEVECTOR_58(myMPC_l, 1);
myMPC_LA_INITIALIZEVECTOR_58(myMPC_s, 1);
info->mu = 0;
myMPC_LA_DOTACC_58(myMPC_l, myMPC_s, &info->mu);
info->mu /= 58;
PRINTTEXT("This is myMPC, a solver generated by FORCES.\n");
PRINTTEXT("(c) Alexander Domahidi, Automatic Control Laboratory, ETH Zurich, 2011-2012.\n");
PRINTTEXT("\n  #it  res_eq   res_ineq     pobj         dobj       dgap     rdgap     mu\n");
PRINTTEXT("  ---------------------------------------------------------------------------\n");
while( 1 ){
info->pobj = 0;
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f0, myMPC_z0, myMPC_grad_cost0, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f1, myMPC_z1, myMPC_grad_cost1, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f2, myMPC_z2, myMPC_grad_cost2, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f3, myMPC_z3, myMPC_grad_cost3, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f4, myMPC_z4, myMPC_grad_cost4, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f5, myMPC_z5, myMPC_grad_cost5, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f6, myMPC_z6, myMPC_grad_cost6, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f7, myMPC_z7, myMPC_grad_cost7, &info->pobj);
myMPC_LA_DENSE_QUADFCN_3(myMPC_H0, myMPC_f8, myMPC_z8, myMPC_grad_cost8, &info->pobj);
myMPC_LA_DENSE_QUADFCN_2(myMPC_H9, myMPC_f9, myMPC_z9, myMPC_grad_cost9, &info->pobj);
info->res_eq = 0;
info->dgap = 0;
myMPC_LA_DENSE_2MVMSUB_4_3_3(myMPC_C0, myMPC_z0, myMPC_D1, myMPC_z1, params->z1, myMPC_v0, myMPC_re0, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z1, myMPC_D2, myMPC_z2, myMPC_c1, myMPC_v1, myMPC_re1, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z2, myMPC_D2, myMPC_z3, myMPC_c2, myMPC_v2, myMPC_re2, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z3, myMPC_D2, myMPC_z4, myMPC_c3, myMPC_v3, myMPC_re3, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z4, myMPC_D2, myMPC_z5, myMPC_c4, myMPC_v4, myMPC_re4, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z5, myMPC_D2, myMPC_z6, myMPC_c5, myMPC_v5, myMPC_re5, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z6, myMPC_D2, myMPC_z7, myMPC_c6, myMPC_v6, myMPC_re6, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_3(myMPC_C1, myMPC_z7, myMPC_D2, myMPC_z8, myMPC_c7, myMPC_v7, myMPC_re7, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_2MVMSUB_2_3_2(myMPC_C1, myMPC_z8, myMPC_D9, myMPC_z9, myMPC_c8, myMPC_v8, myMPC_re8, &info->dgap, &info->res_eq);
myMPC_LA_DENSE_MTVM_4_3(myMPC_C0, myMPC_v0, myMPC_grad_eq0);
myMPC_LA_DENSE_MTVM2_2_3_4(myMPC_C1, myMPC_v1, myMPC_D1, myMPC_v0, myMPC_grad_eq1);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v2, myMPC_D2, myMPC_v1, myMPC_grad_eq2);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v3, myMPC_D2, myMPC_v2, myMPC_grad_eq3);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v4, myMPC_D2, myMPC_v3, myMPC_grad_eq4);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v5, myMPC_D2, myMPC_v4, myMPC_grad_eq5);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v6, myMPC_D2, myMPC_v5, myMPC_grad_eq6);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v7, myMPC_D2, myMPC_v6, myMPC_grad_eq7);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_v8, myMPC_D2, myMPC_v7, myMPC_grad_eq8);
myMPC_LA_DENSE_MTVM_2_2(myMPC_D9, myMPC_v8, myMPC_grad_eq9);
info->res_ineq = 0;
myMPC_LA_VSUBADD3_3(myMPC_lb0, myMPC_z0, myMPC_lbIdx0, myMPC_llb0, myMPC_slb0, myMPC_rilb0, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z0, myMPC_ubIdx0, myMPC_ub0, myMPC_lub0, myMPC_sub0, myMPC_riub0, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb1, myMPC_z1, myMPC_lbIdx1, myMPC_llb1, myMPC_slb1, myMPC_rilb1, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z1, myMPC_ubIdx1, myMPC_ub1, myMPC_lub1, myMPC_sub1, myMPC_riub1, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb2, myMPC_z2, myMPC_lbIdx2, myMPC_llb2, myMPC_slb2, myMPC_rilb2, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z2, myMPC_ubIdx2, myMPC_ub2, myMPC_lub2, myMPC_sub2, myMPC_riub2, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb3, myMPC_z3, myMPC_lbIdx3, myMPC_llb3, myMPC_slb3, myMPC_rilb3, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z3, myMPC_ubIdx3, myMPC_ub3, myMPC_lub3, myMPC_sub3, myMPC_riub3, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb4, myMPC_z4, myMPC_lbIdx4, myMPC_llb4, myMPC_slb4, myMPC_rilb4, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z4, myMPC_ubIdx4, myMPC_ub4, myMPC_lub4, myMPC_sub4, myMPC_riub4, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb5, myMPC_z5, myMPC_lbIdx5, myMPC_llb5, myMPC_slb5, myMPC_rilb5, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z5, myMPC_ubIdx5, myMPC_ub5, myMPC_lub5, myMPC_sub5, myMPC_riub5, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb6, myMPC_z6, myMPC_lbIdx6, myMPC_llb6, myMPC_slb6, myMPC_rilb6, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z6, myMPC_ubIdx6, myMPC_ub6, myMPC_lub6, myMPC_sub6, myMPC_riub6, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb7, myMPC_z7, myMPC_lbIdx7, myMPC_llb7, myMPC_slb7, myMPC_rilb7, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z7, myMPC_ubIdx7, myMPC_ub7, myMPC_lub7, myMPC_sub7, myMPC_riub7, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_3(myMPC_lb8, myMPC_z8, myMPC_lbIdx8, myMPC_llb8, myMPC_slb8, myMPC_rilb8, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_3(myMPC_z8, myMPC_ubIdx8, myMPC_ub8, myMPC_lub8, myMPC_sub8, myMPC_riub8, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD3_2(myMPC_lb9, myMPC_z9, myMPC_lbIdx9, myMPC_llb9, myMPC_slb9, myMPC_rilb9, &info->dgap, &info->res_ineq);
myMPC_LA_VSUBADD2_2(myMPC_z9, myMPC_ubIdx9, myMPC_ub9, myMPC_lub9, myMPC_sub9, myMPC_riub9, &info->dgap, &info->res_ineq);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub0, myMPC_sub0, myMPC_riub0, myMPC_llb0, myMPC_slb0, myMPC_rilb0, myMPC_lbIdx0, myMPC_ubIdx0, myMPC_grad_ineq0, myMPC_lubbysub0, myMPC_llbbyslb0);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub1, myMPC_sub1, myMPC_riub1, myMPC_llb1, myMPC_slb1, myMPC_rilb1, myMPC_lbIdx1, myMPC_ubIdx1, myMPC_grad_ineq1, myMPC_lubbysub1, myMPC_llbbyslb1);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub2, myMPC_sub2, myMPC_riub2, myMPC_llb2, myMPC_slb2, myMPC_rilb2, myMPC_lbIdx2, myMPC_ubIdx2, myMPC_grad_ineq2, myMPC_lubbysub2, myMPC_llbbyslb2);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub3, myMPC_sub3, myMPC_riub3, myMPC_llb3, myMPC_slb3, myMPC_rilb3, myMPC_lbIdx3, myMPC_ubIdx3, myMPC_grad_ineq3, myMPC_lubbysub3, myMPC_llbbyslb3);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub4, myMPC_sub4, myMPC_riub4, myMPC_llb4, myMPC_slb4, myMPC_rilb4, myMPC_lbIdx4, myMPC_ubIdx4, myMPC_grad_ineq4, myMPC_lubbysub4, myMPC_llbbyslb4);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub5, myMPC_sub5, myMPC_riub5, myMPC_llb5, myMPC_slb5, myMPC_rilb5, myMPC_lbIdx5, myMPC_ubIdx5, myMPC_grad_ineq5, myMPC_lubbysub5, myMPC_llbbyslb5);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub6, myMPC_sub6, myMPC_riub6, myMPC_llb6, myMPC_slb6, myMPC_rilb6, myMPC_lbIdx6, myMPC_ubIdx6, myMPC_grad_ineq6, myMPC_lubbysub6, myMPC_llbbyslb6);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub7, myMPC_sub7, myMPC_riub7, myMPC_llb7, myMPC_slb7, myMPC_rilb7, myMPC_lbIdx7, myMPC_ubIdx7, myMPC_grad_ineq7, myMPC_lubbysub7, myMPC_llbbyslb7);
myMPC_LA_INEQ_B_GRAD_3_3_3(myMPC_lub8, myMPC_sub8, myMPC_riub8, myMPC_llb8, myMPC_slb8, myMPC_rilb8, myMPC_lbIdx8, myMPC_ubIdx8, myMPC_grad_ineq8, myMPC_lubbysub8, myMPC_llbbyslb8);
myMPC_LA_INEQ_B_GRAD_2_2_2(myMPC_lub9, myMPC_sub9, myMPC_riub9, myMPC_llb9, myMPC_slb9, myMPC_rilb9, myMPC_lbIdx9, myMPC_ubIdx9, myMPC_grad_ineq9, myMPC_lubbysub9, myMPC_llbbyslb9);
info->dobj = info->pobj - info->dgap;
info->rdgap = info->pobj ? info->dgap / info->pobj : 1e6;
if( info->rdgap < 0 ) info->rdgap = -info->rdgap;
PRINTTEXT("  %3d  %3.1e  %3.1e  %+6.4e  %+6.4e  %+3.1e  %3.1e  %3.1e\n",info->it, info->res_eq, info->res_ineq, info->pobj, info->dobj, info->dgap, info->rdgap, info->mu);
if( info->mu < myMPC_SET_ACC_KKTCOMPL
    && (info->rdgap < myMPC_SET_ACC_RDGAP || info->dgap < myMPC_SET_ACC_KKTCOMPL)
    && info->res_eq < myMPC_SET_ACC_RESEQ
    && info->res_ineq < myMPC_SET_ACC_RESINEQ ){
PRINTTEXT("OPTIMAL (within RESEQ=%2.1e, RESINEQ=%2.1e, (R)DGAP=(%2.1e)%2.1e, MU=%2.1e).\n",myMPC_SET_ACC_RESEQ, myMPC_SET_ACC_RESINEQ,myMPC_SET_ACC_KKTCOMPL,myMPC_SET_ACC_RDGAP,myMPC_SET_ACC_KKTCOMPL);
exitcode = myMPC_OPTIMAL; break; }
if( info->it == myMPC_SET_MAXIT ){
PRINTTEXT("Maximum number of iterations reached, exiting.\n");
exitcode = myMPC_MAXITREACHED; break; }
myMPC_LA_VVADD3_29(myMPC_grad_cost, myMPC_grad_eq, myMPC_grad_ineq, myMPC_rd);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb0, myMPC_lbIdx0, myMPC_lubbysub0, myMPC_ubIdx0, myMPC_Phi0);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi0);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb1, myMPC_lbIdx1, myMPC_lubbysub1, myMPC_ubIdx1, myMPC_Phi1);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi1);
myMPC_LA_DENSE_MATRIXFORWARDSUB_4_3(myMPC_Phi0, myMPC_C0, myMPC_V0);
myMPC_LA_DENSE_MATRIXFORWARDSUB_4_3(myMPC_Phi1, myMPC_D1, myMPC_W1);
myMPC_LA_DENSE_MMT2_4_3_3(myMPC_V0, myMPC_W1, myMPC_Yd0);
myMPC_LA_DENSE_CHOL_4(myMPC_Yd0, myMPC_Ld0);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi0, myMPC_rd0, myMPC_Lbyrd0);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi1, myMPC_rd1, myMPC_Lbyrd1);
myMPC_LA_DENSE_2MVMSUB2_4_3_3(myMPC_V0, myMPC_Lbyrd0, myMPC_W1, myMPC_Lbyrd1, myMPC_re0, myMPC_beta0);
myMPC_LA_DENSE_FORWARDSUB_4(myMPC_Ld0, myMPC_beta0, myMPC_yy0);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb2, myMPC_lbIdx2, myMPC_lubbysub2, myMPC_ubIdx2, myMPC_Phi2);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi2);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi1, myMPC_C1, myMPC_V1);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi2, myMPC_D2, myMPC_W2);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V1, myMPC_W2, myMPC_Yd1);
myMPC_LA_DENSE_MMTM_4_3_2(myMPC_W1, myMPC_V1, myMPC_Ysd1);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_4(myMPC_Ld0, myMPC_Ysd1, myMPC_Lsd1);
myMPC_LA_DENSE_MMTSUB_2_4(myMPC_Lsd1, myMPC_Yd1);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd1, myMPC_Ld1);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi2, myMPC_rd2, myMPC_Lbyrd2);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi1, myMPC_rd1, myMPC_Lbyrd1);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V1, myMPC_Lbyrd1, myMPC_W2, myMPC_Lbyrd2, myMPC_re1, myMPC_beta1);
myMPC_LA_DENSE_MVMSUB1_2_4(myMPC_Lsd1, myMPC_yy0, myMPC_beta1, myMPC_bmy1);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld1, myMPC_bmy1, myMPC_yy1);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb3, myMPC_lbIdx3, myMPC_lubbysub3, myMPC_ubIdx3, myMPC_Phi3);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi3);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi2, myMPC_C1, myMPC_V2);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi3, myMPC_D2, myMPC_W3);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V2, myMPC_W3, myMPC_Yd2);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W2, myMPC_V2, myMPC_Ysd2);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld1, myMPC_Ysd2, myMPC_Lsd2);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd2, myMPC_Yd2);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd2, myMPC_Ld2);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi3, myMPC_rd3, myMPC_Lbyrd3);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi2, myMPC_rd2, myMPC_Lbyrd2);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V2, myMPC_Lbyrd2, myMPC_W3, myMPC_Lbyrd3, myMPC_re2, myMPC_beta2);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd2, myMPC_yy1, myMPC_beta2, myMPC_bmy2);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld2, myMPC_bmy2, myMPC_yy2);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb4, myMPC_lbIdx4, myMPC_lubbysub4, myMPC_ubIdx4, myMPC_Phi4);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi4);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi3, myMPC_C1, myMPC_V3);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi4, myMPC_D2, myMPC_W4);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V3, myMPC_W4, myMPC_Yd3);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W3, myMPC_V3, myMPC_Ysd3);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld2, myMPC_Ysd3, myMPC_Lsd3);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd3, myMPC_Yd3);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd3, myMPC_Ld3);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi4, myMPC_rd4, myMPC_Lbyrd4);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi3, myMPC_rd3, myMPC_Lbyrd3);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V3, myMPC_Lbyrd3, myMPC_W4, myMPC_Lbyrd4, myMPC_re3, myMPC_beta3);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd3, myMPC_yy2, myMPC_beta3, myMPC_bmy3);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld3, myMPC_bmy3, myMPC_yy3);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb5, myMPC_lbIdx5, myMPC_lubbysub5, myMPC_ubIdx5, myMPC_Phi5);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi5);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi4, myMPC_C1, myMPC_V4);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi5, myMPC_D2, myMPC_W5);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V4, myMPC_W5, myMPC_Yd4);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W4, myMPC_V4, myMPC_Ysd4);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld3, myMPC_Ysd4, myMPC_Lsd4);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd4, myMPC_Yd4);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd4, myMPC_Ld4);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi5, myMPC_rd5, myMPC_Lbyrd5);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi4, myMPC_rd4, myMPC_Lbyrd4);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V4, myMPC_Lbyrd4, myMPC_W5, myMPC_Lbyrd5, myMPC_re4, myMPC_beta4);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd4, myMPC_yy3, myMPC_beta4, myMPC_bmy4);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld4, myMPC_bmy4, myMPC_yy4);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb6, myMPC_lbIdx6, myMPC_lubbysub6, myMPC_ubIdx6, myMPC_Phi6);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi6);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi5, myMPC_C1, myMPC_V5);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi6, myMPC_D2, myMPC_W6);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V5, myMPC_W6, myMPC_Yd5);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W5, myMPC_V5, myMPC_Ysd5);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld4, myMPC_Ysd5, myMPC_Lsd5);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd5, myMPC_Yd5);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd5, myMPC_Ld5);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi6, myMPC_rd6, myMPC_Lbyrd6);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi5, myMPC_rd5, myMPC_Lbyrd5);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V5, myMPC_Lbyrd5, myMPC_W6, myMPC_Lbyrd6, myMPC_re5, myMPC_beta5);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd5, myMPC_yy4, myMPC_beta5, myMPC_bmy5);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld5, myMPC_bmy5, myMPC_yy5);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb7, myMPC_lbIdx7, myMPC_lubbysub7, myMPC_ubIdx7, myMPC_Phi7);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi7);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi6, myMPC_C1, myMPC_V6);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi7, myMPC_D2, myMPC_W7);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V6, myMPC_W7, myMPC_Yd6);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W6, myMPC_V6, myMPC_Ysd6);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld5, myMPC_Ysd6, myMPC_Lsd6);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd6, myMPC_Yd6);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd6, myMPC_Ld6);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi7, myMPC_rd7, myMPC_Lbyrd7);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi6, myMPC_rd6, myMPC_Lbyrd6);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V6, myMPC_Lbyrd6, myMPC_W7, myMPC_Lbyrd7, myMPC_re6, myMPC_beta6);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd6, myMPC_yy5, myMPC_beta6, myMPC_bmy6);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld6, myMPC_bmy6, myMPC_yy6);
myMPC_LA_INEQ_DENSE_HESS_3_3_3(myMPC_H0, myMPC_llbbyslb8, myMPC_lbIdx8, myMPC_lubbysub8, myMPC_ubIdx8, myMPC_Phi8);
myMPC_LA_DENSE_CHOL2_3(myMPC_Phi8);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi7, myMPC_C1, myMPC_V7);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi8, myMPC_D2, myMPC_W8);
myMPC_LA_DENSE_MMT2_2_3_3(myMPC_V7, myMPC_W8, myMPC_Yd7);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W7, myMPC_V7, myMPC_Ysd7);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld6, myMPC_Ysd7, myMPC_Lsd7);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd7, myMPC_Yd7);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd7, myMPC_Ld7);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi8, myMPC_rd8, myMPC_Lbyrd8);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi7, myMPC_rd7, myMPC_Lbyrd7);
myMPC_LA_DENSE_2MVMSUB2_2_3_3(myMPC_V7, myMPC_Lbyrd7, myMPC_W8, myMPC_Lbyrd8, myMPC_re7, myMPC_beta7);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd7, myMPC_yy6, myMPC_beta7, myMPC_bmy7);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld7, myMPC_bmy7, myMPC_yy7);
myMPC_LA_INEQ_DENSE_HESS_2_2_2(myMPC_H9, myMPC_llbbyslb9, myMPC_lbIdx9, myMPC_lubbysub9, myMPC_ubIdx9, myMPC_Phi9);
myMPC_LA_DENSE_CHOL2_2(myMPC_Phi9);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_3(myMPC_Phi8, myMPC_C1, myMPC_V8);
myMPC_LA_DENSE_MATRIXFORWARDSUB_2_2(myMPC_Phi9, myMPC_D9, myMPC_W9);
myMPC_LA_DENSE_MMT2_2_3_2(myMPC_V8, myMPC_W9, myMPC_Yd8);
myMPC_LA_DENSE_MMTM_2_3_2(myMPC_W8, myMPC_V8, myMPC_Ysd8);
myMPC_LA_DENSE_MATRIXTFORWARDSUB_2_2(myMPC_Ld7, myMPC_Ysd8, myMPC_Lsd8);
myMPC_LA_DENSE_MMTSUB_2_2(myMPC_Lsd8, myMPC_Yd8);
myMPC_LA_DENSE_CHOL_2(myMPC_Yd8, myMPC_Ld8);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Phi9, myMPC_rd9, myMPC_Lbyrd9);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi8, myMPC_rd8, myMPC_Lbyrd8);
myMPC_LA_DENSE_2MVMSUB2_2_3_2(myMPC_V8, myMPC_Lbyrd8, myMPC_W9, myMPC_Lbyrd9, myMPC_re8, myMPC_beta8);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd8, myMPC_yy7, myMPC_beta8, myMPC_bmy8);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld8, myMPC_bmy8, myMPC_yy8);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld8, myMPC_yy8, myMPC_dvaff8);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd8, myMPC_dvaff8, myMPC_yy7, myMPC_bmy7);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld7, myMPC_bmy7, myMPC_dvaff7);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd7, myMPC_dvaff7, myMPC_yy6, myMPC_bmy6);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld6, myMPC_bmy6, myMPC_dvaff6);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd6, myMPC_dvaff6, myMPC_yy5, myMPC_bmy5);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld5, myMPC_bmy5, myMPC_dvaff5);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd5, myMPC_dvaff5, myMPC_yy4, myMPC_bmy4);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld4, myMPC_bmy4, myMPC_dvaff4);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd4, myMPC_dvaff4, myMPC_yy3, myMPC_bmy3);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld3, myMPC_bmy3, myMPC_dvaff3);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd3, myMPC_dvaff3, myMPC_yy2, myMPC_bmy2);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld2, myMPC_bmy2, myMPC_dvaff2);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd2, myMPC_dvaff2, myMPC_yy1, myMPC_bmy1);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld1, myMPC_bmy1, myMPC_dvaff1);
myMPC_LA_DENSE_MTVMSUB_2_4(myMPC_Lsd1, myMPC_dvaff1, myMPC_yy0, myMPC_bmy0);
myMPC_LA_DENSE_BACKWARDSUB_4(myMPC_Ld0, myMPC_bmy0, myMPC_dvaff0);
myMPC_LA_DENSE_MTVM_4_3(myMPC_C0, myMPC_dvaff0, myMPC_grad_eq0);
myMPC_LA_DENSE_MTVM2_2_3_4(myMPC_C1, myMPC_dvaff1, myMPC_D1, myMPC_dvaff0, myMPC_grad_eq1);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff2, myMPC_D2, myMPC_dvaff1, myMPC_grad_eq2);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff3, myMPC_D2, myMPC_dvaff2, myMPC_grad_eq3);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff4, myMPC_D2, myMPC_dvaff3, myMPC_grad_eq4);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff5, myMPC_D2, myMPC_dvaff4, myMPC_grad_eq5);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff6, myMPC_D2, myMPC_dvaff5, myMPC_grad_eq6);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff7, myMPC_D2, myMPC_dvaff6, myMPC_grad_eq7);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvaff8, myMPC_D2, myMPC_dvaff7, myMPC_grad_eq8);
myMPC_LA_DENSE_MTVM_2_2(myMPC_D9, myMPC_dvaff8, myMPC_grad_eq9);
myMPC_LA_VSUB2_29(myMPC_rd, myMPC_grad_eq, myMPC_rd);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi0, myMPC_rd0, myMPC_dzaff0);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi1, myMPC_rd1, myMPC_dzaff1);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi2, myMPC_rd2, myMPC_dzaff2);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi3, myMPC_rd3, myMPC_dzaff3);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi4, myMPC_rd4, myMPC_dzaff4);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi5, myMPC_rd5, myMPC_dzaff5);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi6, myMPC_rd6, myMPC_dzaff6);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi7, myMPC_rd7, myMPC_dzaff7);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi8, myMPC_rd8, myMPC_dzaff8);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_2(myMPC_Phi9, myMPC_rd9, myMPC_dzaff9);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff0, myMPC_lbIdx0, myMPC_rilb0, myMPC_dslbaff0);
myMPC_LA_VSUB3_3(myMPC_llbbyslb0, myMPC_dslbaff0, myMPC_llb0, myMPC_dllbaff0);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub0, myMPC_dzaff0, myMPC_ubIdx0, myMPC_dsubaff0);
myMPC_LA_VSUB3_3(myMPC_lubbysub0, myMPC_dsubaff0, myMPC_lub0, myMPC_dlubaff0);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff1, myMPC_lbIdx1, myMPC_rilb1, myMPC_dslbaff1);
myMPC_LA_VSUB3_3(myMPC_llbbyslb1, myMPC_dslbaff1, myMPC_llb1, myMPC_dllbaff1);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub1, myMPC_dzaff1, myMPC_ubIdx1, myMPC_dsubaff1);
myMPC_LA_VSUB3_3(myMPC_lubbysub1, myMPC_dsubaff1, myMPC_lub1, myMPC_dlubaff1);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff2, myMPC_lbIdx2, myMPC_rilb2, myMPC_dslbaff2);
myMPC_LA_VSUB3_3(myMPC_llbbyslb2, myMPC_dslbaff2, myMPC_llb2, myMPC_dllbaff2);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub2, myMPC_dzaff2, myMPC_ubIdx2, myMPC_dsubaff2);
myMPC_LA_VSUB3_3(myMPC_lubbysub2, myMPC_dsubaff2, myMPC_lub2, myMPC_dlubaff2);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff3, myMPC_lbIdx3, myMPC_rilb3, myMPC_dslbaff3);
myMPC_LA_VSUB3_3(myMPC_llbbyslb3, myMPC_dslbaff3, myMPC_llb3, myMPC_dllbaff3);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub3, myMPC_dzaff3, myMPC_ubIdx3, myMPC_dsubaff3);
myMPC_LA_VSUB3_3(myMPC_lubbysub3, myMPC_dsubaff3, myMPC_lub3, myMPC_dlubaff3);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff4, myMPC_lbIdx4, myMPC_rilb4, myMPC_dslbaff4);
myMPC_LA_VSUB3_3(myMPC_llbbyslb4, myMPC_dslbaff4, myMPC_llb4, myMPC_dllbaff4);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub4, myMPC_dzaff4, myMPC_ubIdx4, myMPC_dsubaff4);
myMPC_LA_VSUB3_3(myMPC_lubbysub4, myMPC_dsubaff4, myMPC_lub4, myMPC_dlubaff4);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff5, myMPC_lbIdx5, myMPC_rilb5, myMPC_dslbaff5);
myMPC_LA_VSUB3_3(myMPC_llbbyslb5, myMPC_dslbaff5, myMPC_llb5, myMPC_dllbaff5);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub5, myMPC_dzaff5, myMPC_ubIdx5, myMPC_dsubaff5);
myMPC_LA_VSUB3_3(myMPC_lubbysub5, myMPC_dsubaff5, myMPC_lub5, myMPC_dlubaff5);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff6, myMPC_lbIdx6, myMPC_rilb6, myMPC_dslbaff6);
myMPC_LA_VSUB3_3(myMPC_llbbyslb6, myMPC_dslbaff6, myMPC_llb6, myMPC_dllbaff6);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub6, myMPC_dzaff6, myMPC_ubIdx6, myMPC_dsubaff6);
myMPC_LA_VSUB3_3(myMPC_lubbysub6, myMPC_dsubaff6, myMPC_lub6, myMPC_dlubaff6);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff7, myMPC_lbIdx7, myMPC_rilb7, myMPC_dslbaff7);
myMPC_LA_VSUB3_3(myMPC_llbbyslb7, myMPC_dslbaff7, myMPC_llb7, myMPC_dllbaff7);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub7, myMPC_dzaff7, myMPC_ubIdx7, myMPC_dsubaff7);
myMPC_LA_VSUB3_3(myMPC_lubbysub7, myMPC_dsubaff7, myMPC_lub7, myMPC_dlubaff7);
myMPC_LA_VSUB_INDEXED_3(myMPC_dzaff8, myMPC_lbIdx8, myMPC_rilb8, myMPC_dslbaff8);
myMPC_LA_VSUB3_3(myMPC_llbbyslb8, myMPC_dslbaff8, myMPC_llb8, myMPC_dllbaff8);
myMPC_LA_VSUB2_INDEXED_3(myMPC_riub8, myMPC_dzaff8, myMPC_ubIdx8, myMPC_dsubaff8);
myMPC_LA_VSUB3_3(myMPC_lubbysub8, myMPC_dsubaff8, myMPC_lub8, myMPC_dlubaff8);
myMPC_LA_VSUB_INDEXED_2(myMPC_dzaff9, myMPC_lbIdx9, myMPC_rilb9, myMPC_dslbaff9);
myMPC_LA_VSUB3_2(myMPC_llbbyslb9, myMPC_dslbaff9, myMPC_llb9, myMPC_dllbaff9);
myMPC_LA_VSUB2_INDEXED_2(myMPC_riub9, myMPC_dzaff9, myMPC_ubIdx9, myMPC_dsubaff9);
myMPC_LA_VSUB3_2(myMPC_lubbysub9, myMPC_dsubaff9, myMPC_lub9, myMPC_dlubaff9);
info->lsit_aff = myMPC_LINESEARCH_BACKTRACKING_AFFINE(myMPC_l, myMPC_s, myMPC_dl_aff, myMPC_ds_aff, &info->step_aff, &info->mu_aff);
if( info->lsit_aff == myMPC_NOPROGRESS ){
PRINTTEXT("Affine line search could not proceed at iteration %d, exiting.\n",info->it+1);
return myMPC_NOPROGRESS;
}
sigma_3rdroot = info->mu_aff / info->mu;
info->sigma = sigma_3rdroot*sigma_3rdroot*sigma_3rdroot;
musigma = info->mu * info->sigma;
myMPC_LA_VSUB5_58(myMPC_ds_aff, myMPC_dl_aff, musigma, myMPC_ccrhs);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub0, myMPC_sub0, myMPC_ubIdx0, myMPC_ccrhsl0, myMPC_slb0, myMPC_lbIdx0, myMPC_rd0);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub1, myMPC_sub1, myMPC_ubIdx1, myMPC_ccrhsl1, myMPC_slb1, myMPC_lbIdx1, myMPC_rd1);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi0, myMPC_rd0, myMPC_Lbyrd0);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi1, myMPC_rd1, myMPC_Lbyrd1);
myMPC_LA_DENSE_2MVMADD_4_3_3(myMPC_V0, myMPC_Lbyrd0, myMPC_W1, myMPC_Lbyrd1, myMPC_beta0);
myMPC_LA_DENSE_FORWARDSUB_4(myMPC_Ld0, myMPC_beta0, myMPC_yy0);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub2, myMPC_sub2, myMPC_ubIdx2, myMPC_ccrhsl2, myMPC_slb2, myMPC_lbIdx2, myMPC_rd2);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi2, myMPC_rd2, myMPC_Lbyrd2);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi1, myMPC_rd1, myMPC_Lbyrd1);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V1, myMPC_Lbyrd1, myMPC_W2, myMPC_Lbyrd2, myMPC_beta1);
myMPC_LA_DENSE_MVMSUB1_2_4(myMPC_Lsd1, myMPC_yy0, myMPC_beta1, myMPC_bmy1);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld1, myMPC_bmy1, myMPC_yy1);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub3, myMPC_sub3, myMPC_ubIdx3, myMPC_ccrhsl3, myMPC_slb3, myMPC_lbIdx3, myMPC_rd3);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi3, myMPC_rd3, myMPC_Lbyrd3);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi2, myMPC_rd2, myMPC_Lbyrd2);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V2, myMPC_Lbyrd2, myMPC_W3, myMPC_Lbyrd3, myMPC_beta2);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd2, myMPC_yy1, myMPC_beta2, myMPC_bmy2);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld2, myMPC_bmy2, myMPC_yy2);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub4, myMPC_sub4, myMPC_ubIdx4, myMPC_ccrhsl4, myMPC_slb4, myMPC_lbIdx4, myMPC_rd4);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi4, myMPC_rd4, myMPC_Lbyrd4);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi3, myMPC_rd3, myMPC_Lbyrd3);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V3, myMPC_Lbyrd3, myMPC_W4, myMPC_Lbyrd4, myMPC_beta3);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd3, myMPC_yy2, myMPC_beta3, myMPC_bmy3);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld3, myMPC_bmy3, myMPC_yy3);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub5, myMPC_sub5, myMPC_ubIdx5, myMPC_ccrhsl5, myMPC_slb5, myMPC_lbIdx5, myMPC_rd5);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi5, myMPC_rd5, myMPC_Lbyrd5);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi4, myMPC_rd4, myMPC_Lbyrd4);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V4, myMPC_Lbyrd4, myMPC_W5, myMPC_Lbyrd5, myMPC_beta4);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd4, myMPC_yy3, myMPC_beta4, myMPC_bmy4);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld4, myMPC_bmy4, myMPC_yy4);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub6, myMPC_sub6, myMPC_ubIdx6, myMPC_ccrhsl6, myMPC_slb6, myMPC_lbIdx6, myMPC_rd6);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi6, myMPC_rd6, myMPC_Lbyrd6);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi5, myMPC_rd5, myMPC_Lbyrd5);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V5, myMPC_Lbyrd5, myMPC_W6, myMPC_Lbyrd6, myMPC_beta5);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd5, myMPC_yy4, myMPC_beta5, myMPC_bmy5);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld5, myMPC_bmy5, myMPC_yy5);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub7, myMPC_sub7, myMPC_ubIdx7, myMPC_ccrhsl7, myMPC_slb7, myMPC_lbIdx7, myMPC_rd7);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi7, myMPC_rd7, myMPC_Lbyrd7);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi6, myMPC_rd6, myMPC_Lbyrd6);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V6, myMPC_Lbyrd6, myMPC_W7, myMPC_Lbyrd7, myMPC_beta6);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd6, myMPC_yy5, myMPC_beta6, myMPC_bmy6);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld6, myMPC_bmy6, myMPC_yy6);
myMPC_LA_VSUB6_INDEXED_3_3_3(myMPC_ccrhsub8, myMPC_sub8, myMPC_ubIdx8, myMPC_ccrhsl8, myMPC_slb8, myMPC_lbIdx8, myMPC_rd8);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi8, myMPC_rd8, myMPC_Lbyrd8);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi7, myMPC_rd7, myMPC_Lbyrd7);
myMPC_LA_DENSE_2MVMADD_2_3_3(myMPC_V7, myMPC_Lbyrd7, myMPC_W8, myMPC_Lbyrd8, myMPC_beta7);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd7, myMPC_yy6, myMPC_beta7, myMPC_bmy7);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld7, myMPC_bmy7, myMPC_yy7);
myMPC_LA_VSUB6_INDEXED_2_2_2(myMPC_ccrhsub9, myMPC_sub9, myMPC_ubIdx9, myMPC_ccrhsl9, myMPC_slb9, myMPC_lbIdx9, myMPC_rd9);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Phi9, myMPC_rd9, myMPC_Lbyrd9);
myMPC_LA_DENSE_FORWARDSUB_3(myMPC_Phi8, myMPC_rd8, myMPC_Lbyrd8);
myMPC_LA_DENSE_2MVMADD_2_3_2(myMPC_V8, myMPC_Lbyrd8, myMPC_W9, myMPC_Lbyrd9, myMPC_beta8);
myMPC_LA_DENSE_MVMSUB1_2_2(myMPC_Lsd8, myMPC_yy7, myMPC_beta8, myMPC_bmy8);
myMPC_LA_DENSE_FORWARDSUB_2(myMPC_Ld8, myMPC_bmy8, myMPC_yy8);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld8, myMPC_yy8, myMPC_dvcc8);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd8, myMPC_dvcc8, myMPC_yy7, myMPC_bmy7);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld7, myMPC_bmy7, myMPC_dvcc7);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd7, myMPC_dvcc7, myMPC_yy6, myMPC_bmy6);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld6, myMPC_bmy6, myMPC_dvcc6);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd6, myMPC_dvcc6, myMPC_yy5, myMPC_bmy5);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld5, myMPC_bmy5, myMPC_dvcc5);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd5, myMPC_dvcc5, myMPC_yy4, myMPC_bmy4);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld4, myMPC_bmy4, myMPC_dvcc4);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd4, myMPC_dvcc4, myMPC_yy3, myMPC_bmy3);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld3, myMPC_bmy3, myMPC_dvcc3);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd3, myMPC_dvcc3, myMPC_yy2, myMPC_bmy2);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld2, myMPC_bmy2, myMPC_dvcc2);
myMPC_LA_DENSE_MTVMSUB_2_2(myMPC_Lsd2, myMPC_dvcc2, myMPC_yy1, myMPC_bmy1);
myMPC_LA_DENSE_BACKWARDSUB_2(myMPC_Ld1, myMPC_bmy1, myMPC_dvcc1);
myMPC_LA_DENSE_MTVMSUB_2_4(myMPC_Lsd1, myMPC_dvcc1, myMPC_yy0, myMPC_bmy0);
myMPC_LA_DENSE_BACKWARDSUB_4(myMPC_Ld0, myMPC_bmy0, myMPC_dvcc0);
myMPC_LA_DENSE_MTVM_4_3(myMPC_C0, myMPC_dvcc0, myMPC_grad_eq0);
myMPC_LA_DENSE_MTVM2_2_3_4(myMPC_C1, myMPC_dvcc1, myMPC_D1, myMPC_dvcc0, myMPC_grad_eq1);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc2, myMPC_D2, myMPC_dvcc1, myMPC_grad_eq2);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc3, myMPC_D2, myMPC_dvcc2, myMPC_grad_eq3);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc4, myMPC_D2, myMPC_dvcc3, myMPC_grad_eq4);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc5, myMPC_D2, myMPC_dvcc4, myMPC_grad_eq5);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc6, myMPC_D2, myMPC_dvcc5, myMPC_grad_eq6);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc7, myMPC_D2, myMPC_dvcc6, myMPC_grad_eq7);
myMPC_LA_DENSE_MTVM2_2_3_2(myMPC_C1, myMPC_dvcc8, myMPC_D2, myMPC_dvcc7, myMPC_grad_eq8);
myMPC_LA_DENSE_MTVM_2_2(myMPC_D9, myMPC_dvcc8, myMPC_grad_eq9);
myMPC_LA_VSUB_29(myMPC_rd, myMPC_grad_eq, myMPC_rd);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi0, myMPC_rd0, myMPC_dzcc0);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi1, myMPC_rd1, myMPC_dzcc1);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi2, myMPC_rd2, myMPC_dzcc2);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi3, myMPC_rd3, myMPC_dzcc3);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi4, myMPC_rd4, myMPC_dzcc4);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi5, myMPC_rd5, myMPC_dzcc5);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi6, myMPC_rd6, myMPC_dzcc6);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi7, myMPC_rd7, myMPC_dzcc7);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_3(myMPC_Phi8, myMPC_rd8, myMPC_dzcc8);
myMPC_LA_DENSE_FORWARDBACKWARDSUB_2(myMPC_Phi9, myMPC_rd9, myMPC_dzcc9);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl0, myMPC_slb0, myMPC_llbbyslb0, myMPC_dzcc0, myMPC_lbIdx0, myMPC_dllbcc0);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub0, myMPC_sub0, myMPC_lubbysub0, myMPC_dzcc0, myMPC_ubIdx0, myMPC_dlubcc0);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl1, myMPC_slb1, myMPC_llbbyslb1, myMPC_dzcc1, myMPC_lbIdx1, myMPC_dllbcc1);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub1, myMPC_sub1, myMPC_lubbysub1, myMPC_dzcc1, myMPC_ubIdx1, myMPC_dlubcc1);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl2, myMPC_slb2, myMPC_llbbyslb2, myMPC_dzcc2, myMPC_lbIdx2, myMPC_dllbcc2);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub2, myMPC_sub2, myMPC_lubbysub2, myMPC_dzcc2, myMPC_ubIdx2, myMPC_dlubcc2);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl3, myMPC_slb3, myMPC_llbbyslb3, myMPC_dzcc3, myMPC_lbIdx3, myMPC_dllbcc3);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub3, myMPC_sub3, myMPC_lubbysub3, myMPC_dzcc3, myMPC_ubIdx3, myMPC_dlubcc3);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl4, myMPC_slb4, myMPC_llbbyslb4, myMPC_dzcc4, myMPC_lbIdx4, myMPC_dllbcc4);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub4, myMPC_sub4, myMPC_lubbysub4, myMPC_dzcc4, myMPC_ubIdx4, myMPC_dlubcc4);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl5, myMPC_slb5, myMPC_llbbyslb5, myMPC_dzcc5, myMPC_lbIdx5, myMPC_dllbcc5);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub5, myMPC_sub5, myMPC_lubbysub5, myMPC_dzcc5, myMPC_ubIdx5, myMPC_dlubcc5);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl6, myMPC_slb6, myMPC_llbbyslb6, myMPC_dzcc6, myMPC_lbIdx6, myMPC_dllbcc6);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub6, myMPC_sub6, myMPC_lubbysub6, myMPC_dzcc6, myMPC_ubIdx6, myMPC_dlubcc6);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl7, myMPC_slb7, myMPC_llbbyslb7, myMPC_dzcc7, myMPC_lbIdx7, myMPC_dllbcc7);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub7, myMPC_sub7, myMPC_lubbysub7, myMPC_dzcc7, myMPC_ubIdx7, myMPC_dlubcc7);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_3(myMPC_ccrhsl8, myMPC_slb8, myMPC_llbbyslb8, myMPC_dzcc8, myMPC_lbIdx8, myMPC_dllbcc8);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_3(myMPC_ccrhsub8, myMPC_sub8, myMPC_lubbysub8, myMPC_dzcc8, myMPC_ubIdx8, myMPC_dlubcc8);
myMPC_LA_VEC_DIVSUB_MULTSUB_INDEXED_2(myMPC_ccrhsl9, myMPC_slb9, myMPC_llbbyslb9, myMPC_dzcc9, myMPC_lbIdx9, myMPC_dllbcc9);
myMPC_LA_VEC_DIVSUB_MULTADD_INDEXED_2(myMPC_ccrhsub9, myMPC_sub9, myMPC_lubbysub9, myMPC_dzcc9, myMPC_ubIdx9, myMPC_dlubcc9);
myMPC_LA_VSUB7_58(myMPC_l, myMPC_ccrhs, myMPC_s, myMPC_dl_cc, myMPC_ds_cc);
myMPC_LA_VADD_29(myMPC_dz_cc, myMPC_dz_aff);
myMPC_LA_VADD_20(myMPC_dv_cc, myMPC_dv_aff);
myMPC_LA_VADD_58(myMPC_dl_cc, myMPC_dl_aff);
myMPC_LA_VADD_58(myMPC_ds_cc, myMPC_ds_aff);
info->lsit_cc = myMPC_LINESEARCH_BACKTRACKING_COMBINED(myMPC_z, myMPC_v, myMPC_l, myMPC_s, myMPC_dz_cc, myMPC_dv_cc, myMPC_dl_cc, myMPC_ds_cc, &info->step_cc, &info->mu);
if( info->lsit_cc == myMPC_NOPROGRESS ){
PRINTTEXT("Line search could not proceed at iteration %d, exiting.\n",info->it+1);
exitcode = myMPC_NOPROGRESS; break;
}
info->it++;
}
output->u1[0] = myMPC_z0[2];

#if myMPC_SET_PRINTLEVEL > 0
info->solvetime = myMPC_toc(&solvertimer);
if( info->it > 1 ){
	PRINTTEXT("Solve time: %5.3f ms (%d iterations)\n\n", info->it, info->solvetime*1000);
} else {
	PRINTTEXT("Solve time: %5.3f ms (%d iteration)\n\n", info->it, info->solvetime*1000);
}
#else
info->solvetime = -1;
#endif
return exitcode;
}
