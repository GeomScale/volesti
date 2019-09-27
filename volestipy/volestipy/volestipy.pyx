#!python
#cython: language_level=3
#cython: boundscheck=False
#cython: wraparound=False

import os
import sys


import numpy as np
cimport numpy as np

def get_time_seed():
   import random
   import time
   return int(time.time())

cdef extern from "bindings.h":
   cdef cppclass HPolytopeCPP:
      HPolytopeCPP() except +
      HPolytopeCPP(double  *A, double *b, int n_varibles, int n_hyperplanes) except +
      double compute_volume(int walk_len, double epsilon, np.npy_int32 seed)
      double generate_samples(int walk_len, int n_samples, double* samples, np.npy_int32 seed)


cdef class HPolytope:
   cdef HPolytopeCPP polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b

   def __cinit__(self, double[:,::1] A, double[::1] b):
      self._A = A
      self._b = b
      # assert n_hyperplanes==b.shape[0]
      # Do we also need to pass strdie to make it more error proof?
      # https://stackoverflow.com/questions/34592961/passing-multidimensional-memoryviews-to-c-function
      n_hyperplanes, n_variables = A.shape[0], A.shape[1]
      self.polytope_cpp = HPolytopeCPP(&A[0,0], &b[0], n_hyperplanes, n_variables)

   def compute_volume(self, int walk_len=2, double epsilon=0.05, method="gaussian_annealing", np.npy_int32 seed=get_time_seed()):
      if method=="gaussian_annealing":
         return self.polytope_cpp.compute_volume(walk_len, epsilon, seed)
      else:
         raise Exception('"{}" is not implemented to compute volume'.format(method))

   def generate_samples(self, int walk_len=2, int n_samples=1000, np.npy_int32 seed=get_time_seed()):
         n_variables = self._A.shape[1]
         cdef double[:,::1] samples = np.zeros((n_samples,  n_variables), dtype=np.float64, order="C")
         self.polytope_cpp.generate_samples(walk_len, n_samples, &samples[0,0], seed)
         return np.asarray(samples)



   @property
   def A(self):
      return np.asarray(self._A)
   @property
   def b(self):
      return np.asarray(self._b)
   @property
   def dimensions(self):
      return self._A.shape[1]
