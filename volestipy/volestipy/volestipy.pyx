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

def get_time_seed():
   import random
   import time
   return int(time.time())


cdef extern from "bindings.h":
   
   cdef cppclass HPolytopeCPP:
   
      HPolytopeCPP() except +
      HPolytopeCPP(double *A, double *b, int n_hyperplanes, int n_variables) except +

      double compute_volume(char* vol_method, char* walk_method, int walk_len, double epsilon, int seed);     

      double generate_samples(int walk_len, int number_of_points, int number_of_points_to_burn, bool boundary, \
         bool cdhr, bool rdhr, bool gaussian, bool set_L, bool billiard, bool ball_walk, double a, double L,  double* samples);


# with respect to the "compute_volume" def
volume_methods = ["sequence_of_balls".encode("UTF-8"), "cooling_gaussian".encode("UTF-8"), "cooling_balls".encode("UTF-8")]
walk_methods = ["uniform_ball".encode("UTF-8"), "CDHR".encode("UTF-8"), "RDHR".encode("UTF-8"), "gaussian_ball".encode("UTF-8"), \
                "gaussian_CDHR".encode("UTF-8"), "gaussian_RDHR".encode("UTF-8"), "uniform_ball".encode("UTF-8"), "billiard".encode("UTF-8")]


cdef class HPolytope:
   cdef HPolytopeCPP polytope_cpp
   cdef double[:,::1] _A
   cdef double[::1] _b
   
   def __cinit__(self, double[:,::1] A, double[::1] b):
      self._A = A
      self._b = b
      n_hyperplanes, n_variables = A.shape[0], A.shape[1]
      self.polytope_cpp = HPolytopeCPP(&A[0,0], &b[0], n_hyperplanes, n_variables)

   def compute_volume(self, walk_len=2, epsilon=0.05, vol_method="sequence_of_balls", walk_method="uniform_ball", \
      np.npy_int32 seed=get_time_seed()):
      vol_method = vol_method.encode("UTF-8")
      walk_method = walk_method.encode("UTF-8")
      if vol_method  in volume_methods:
         if walk_method in walk_methods:
            return self.polytope_cpp.compute_volume(vol_method, walk_method, walk_len, epsilon, seed)
         else:
            raise Exception('"{}" is not implemented to walk methods. Available methods are: {}'.format(walk_method, walk_methods))
      else:
         raise Exception('"{}" is not implemented to compute volume. Available methods are: {}'.format(vol_method, volume_methods))


# we need to build a Python function for getting a starting point depending on the polytope 
   def generate_samples(self, walk_len=1, number_of_points=1000, number_of_points_to_burn=100, boundary=False, cdhr=True, rdhr=False, gaussian=False, set_L=False, billiard=False, ball_walk=False, a=0, L=0):
      n_variables = self._A.shape[1]
      cdef double[:,::1] samples = np.zeros((number_of_points,  n_variables), dtype=np.float64, order="C")
      
      print("volestipy.pyx: This is just before I call for the generate_samples function \n\n")      
      
      self.polytope_cpp.generate_samples(walk_len, number_of_points, number_of_points_to_burn, boundary, cdhr, rdhr, gaussian, set_L, billiard, ball_walk, a, L, &samples[0,0])
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
