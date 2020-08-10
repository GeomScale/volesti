#!python
#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False

import os
import sys
import numpy as np
cimport numpy as np
from libcpp cimport bool
from cpython cimport bool

# set the time
def get_time_seed():
   import random
   import time
   return int(time.time())

# get classes from the bindings.h file
cdef extern from "bindings.h":
   
   # the HPolytopeCPP class along with its functions
   cdef cppclass HPolytopeCPP:
    
      # initialization 
      HPolytopeCPP() except +
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables) except +

      # compute volume
      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);

      # random sampling
      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, \
         bool cdhr, bool rdhr, bool gaussian, bool set_L, bool billiard, bool ball_walk, double a, double L,  double* samples);

      # rounding H-Polytope
      void rounding(char* rounding_method, double* new_A, double* new_b, double* T_matrix, double* shift, double &round_value);
      
   # the lowDimPolytopeCPP class along with its functions
   cdef cppclass lowDimHPolytopeCPP:
      
      # initialization       
      lowDimHPolytopeCPP() except +
      lowDimHPolytopeCPP(double *A, double *b, double *Aeq, double *beq, int n_rows_of_A, int n_cols_of_A, int n_row_of_Aeq, int n_cols_of_Aeq) except +
      
      # get full dimensional polytope
      int full_dimensiolal_polytope(double* N_extra_trans, double* shift, double* A_full_extra_trans, double* b_full)
      
      

# lists with the methods supported by volesti for volume approximation and random walk
volume_methods = ["sequence_of_balls".encode("UTF-8"), "cooling_gaussian".encode("UTF-8"), "cooling_balls".encode("UTF-8")]
walk_methods = ["uniform_ball".encode("UTF-8"), "CDHR".encode("UTF-8"), "RDHR".encode("UTF-8"), "gaussian_ball".encode("UTF-8"), \
                "gaussian_CDHR".encode("UTF-8"), "gaussian_RDHR".encode("UTF-8"), "uniform_ball".encode("UTF-8"), "billiard".encode("UTF-8")]
rounding_methods = ["min_ellipsoid".encode("UTF-8"), "svd".encode("UTF-8"), "max_ellipsoid".encode("UTF-8")]


# build the HPolytope class - the 'polytope_cpp' is an instance of the HPolytopeCPP class described on the 'bindings.cpp' file
cdef class HPolytope:
   
   cdef HPolytopeCPP polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b

# set the specs of the class
   def __cinit__(self, double[:,::1] A, double[::1] b):
      self._A = A
      self._b = b
      n_hyperplanes, n_variables = A.shape[0], A.shape[1]
      self.polytope_cpp = HPolytopeCPP(&A[0,0], &b[0], n_hyperplanes, n_variables)

#  this is where the volesti functions are getting their python interface; first the compute_volume() function
   def compute_volume(self, walk_len = 2, epsilon = 0.05, vol_method = "sequence_of_balls", walk_method = "uniform_ball", \
      np.npy_int32 seed=get_time_seed()):

      vol_method = vol_method.encode("UTF-8")
      walk_method = walk_method.encode("UTF-8")

      if vol_method in volume_methods:
         if walk_method in walk_methods:
            return self.polytope_cpp.compute_volume(vol_method, walk_method, walk_len, epsilon, seed)
         else:
            raise Exception('"{}" is not implemented to walk methods. Available methods are: {}'.format(walk_method, walk_methods))
      else:
         raise Exception('"{}" is not implemented to compute volume. Available methods are: {}'.format(vol_method, volume_methods))

# likewise, the generate_samples() function
   def generate_samples(self, walk_len = 1, number_of_points = 1000, number_of_points_to_burn = 100, boundary = False, cdhr=True, \
      rdhr = False, gaussian = False, set_L = False, billiard = False, ball_walk = False, a = 0, L = 0):

      n_variables = self._A.shape[1]
      cdef double[:,::1] samples = np.zeros((number_of_points,  n_variables), dtype = np.float64, order = "C")

      self.polytope_cpp.generate_samples(walk_len, number_of_points, number_of_points_to_burn, boundary, cdhr, rdhr, gaussian, set_L, billiard, ball_walk, a, L, &samples[0,0])
      return np.asarray(samples)      # we need to build a Python function for getting a starting point depending on the polytope

# the rounding() function; like the compute_volume; there are more than one methods for this step
   def rounding(self, rounding_method = 'max_ellipsoid'):

      # get the dimensions of the items about to build
      n_hyperplanes, n_variables = self._A.shape[0], self._A.shape[1]

      # set the variables of those items; notice that they are all cdef type except of the last one which is about to be used both as a C++ and a Python variable
      cdef double[:,::1] new_A = np.zeros((n_hyperplanes, n_variables), dtype=np.float64, order="C")
      cdef double[::1] new_b = np.zeros(n_hyperplanes, dtype=np.float64, order="C")
      cdef double[:,::1] T_matrix = np.zeros((n_variables, n_variables), dtype=np.float64, order="C")
      cdef double[::1] shift = np.zeros((n_variables), dtype=np.float64, order="C")
      cdef double round_value

      # transform the rounding_method variable to UTF-8 coding
      rounding_method = rounding_method.encode("UTF-8")

      # check whether the rounding method the user asked for, is actually among those volestipy supports
      if rounding_method in rounding_methods:
         self.polytope_cpp.rounding(rounding_method, &new_A[0,0], &new_b[0], &T_matrix[0,0], &shift[0], round_value)
         return np.asarray(new_A),np.asarray(new_b),np.asarray(T_matrix),np.asarray(shift),np.asarray(round_value)
      else:
         raise Exception('"{}" is not implemented to walk types. Available methods are: {}'.format(rounding_method, rounding_methods))

   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def dimensions(self):
      return self._A.shape[1]



cdef class low_dim_HPolytope:

   cdef lowDimHPolytopeCPP low_dim_polytope_cpp

   cdef double[:,::1] _A
   cdef double[::1] _b
   cdef double [:,::1] _Aeq
   cdef double[::1] _beq

# set the specs of the class
   def __cinit__(self, double[:,::1] A, double[::1] b, double[:,::1] Aeq, double[::1] beq):
      self._A = A
      self._b = b
      self._Aeq = Aeq
      self._beq = beq
      n_rows_of_A, n_cols_of_A = A.shape[0], A.shape[1] 
      n_row_of_Aeq, n_cols_of_Aeq = Aeq.shape[0], Aeq.shape[1] 
      
      self.low_dim_polytope_cpp = lowDimHPolytopeCPP(&A[0,0], &b[0], &Aeq[0,0], &beq[0], n_rows_of_A, n_cols_of_A, n_row_of_Aeq, n_cols_of_Aeq)

   
   # the get_full_dimensional_polytope() function(); that needs to run in case the user does not provide volestipy with a full dimensional polytope
   def full_dimensiolal_polytope(self):
      
      # get dimensions of the initial S (Aeq) matrix
      m = self._Aeq.shape[0]
      n = self._Aeq.shape[1]
      k = self._A.shape[0]
      
      # regarding the actual full_dimensional_polytope() output 
      cdef double[:,::1] N_extra_trans = np.zeros((n, n), dtype=np.float64, order="C")   # the number of lines here are at least n-m; but we do not know their exact number
      cdef double[::1] shift = np.zeros((n), dtype=np.float64, order="C")
      
      # regarding the data needed to build an HPolytopeCPP object with the A and b of the full dimensional polytope
      # that we do not know is the exact number of cols of the A_full and the N matrices, which is the same. As we get their transpose matrices from the .cpp code, that's why we see the unknown in the rows.
      cdef double[:,::1] A_full_extra_trans = np.zeros((n,k), dtype=np.float64, order="C")    # like in the case of the N_trans matrix, we do not know the exact number of lines of this
      cdef double[::1] b_full = np.zeros((k), dtype=np.float64, order="C")
      
      # regarding the actual dimensions of the matrices and vectors of the full dimensional polytope
      cpdef int n_of_cols_in_N
   
      # call the C++ class to get the full_dimensional polytope
      n_of_cols_in_N = self.low_dim_polytope_cpp.full_dimensiolal_polytope(&N_extra_trans[0,0], &shift[0], &A_full_extra_trans[0,0], &b_full[0])
      
      # get a matrix with exactly the number of lines and columns that N expands to
      N = np.zeros((n, n_of_cols_in_N), dtype=np.float64, order="C")
      for i in range(n):
         for j in range(n_of_cols_in_N):
            N[i,j] = np.asarray(N_extra_trans[j,i])
      del N_extra_trans
         
      # and now the full dimensional polytope's specs; A matrix
      # the b vector has not dim issues as we know its length before running the get_full_dimensional_polytope() function
      A_full = np.zeros((k, n_of_cols_in_N), dtype=np.float64, order="C")
      for i in range(k):
         for j in range(n_of_cols_in_N):
            A_full[i,j] = np.asarray(A_full_extra_trans[j,i])
      del A_full_extra_trans
      
      # finally, we need to build an HP object for the full dumensional polytope we got
      full_dimensional_polytope = HPolytope(A_full,b_full)
      
      # and we print the shift and the N matrix and keep the full dim HP object as return
      print("The shift vector is: ")
      print(np.asarray(shift))
      print("------------")

      print("The N matrix is: \n")
      print(np.asarray(N))
      print("-------------")

      return full_dimensional_polytope



   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def Aeq(self):
      return np.asarray(self._Aeq)
   @property
   def beq(self):
      return np.asarray(self._beq)   
   @property
   def dimensions(self):
      return self._A.shape[1]   