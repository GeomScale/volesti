#include "InteriorPointMethod.h"

double t0 = 0.0001;
double nu = 1.5;
double tolerance=0.000001;
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
double IPMprimaldual(mat &A, vec &b, vec &c, vec &p) {
  int max_iter = 150;
  int iter = 0;
  double t = t0;
  double dualf,f;
  vec dualArg;
  while (iter < max_iter) {
    vec temp = c * t;
    NewtonsMethoddual(A, b, temp, p,t,dualf,dualArg,f);
    double dualityGap=f-dualf;
    //cout<<dualityGap<<"\n";
    if(dualityGap<tolerance){
      break;
    }
    t = t * nu;
    iter++;
  }
  return as_scalar(c.t() * p);
}
