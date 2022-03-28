#include "InteriorPointMethod.h"

double t0 = 0.0001;
double nu = 1.1;

double IPM(mat &A, vec &b, vec &c, vec &p) {
  int max_iter = 150;
  int iter = 0;
  double t = t0;
  while (iter < max_iter) {
    vec temp = c * t;
    NewtonsMethod(A, b, temp, p);
    t = t * nu;
    iter++;
  }
  return as_scalar(c.t() * p);
}
