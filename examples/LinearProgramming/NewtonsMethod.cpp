#include "InteriorPointMethod.h"
#include <limits>
double alphaparam=0.1;
double betaparam=0.5;
double threshold = 0.00000001;
double e = numeric_limits<double>::epsilon();
double ComputeF(mat &A, vec &b , vec &c,vec p){
  vec temp=b-A*p;
  temp.transform( [](double val) { return log(max(val,e)); } );

  return as_scalar(c.t()*p-accu(temp));
}

void NewtonsMethod(mat &A, vec &b, vec &c, vec &p) {
  int max_iter = 100;
  int iter = 0;
  int n=A.n_cols;
  int m=A.n_rows;
  vec Grad = vec(n);
  Grad.fill(1);
  mat H = mat(n, n);
  double f;
  while (iter < max_iter && norm(Grad) > threshold) {
    vec temp = b-A * p;
    temp.transform( [](double val) { return 1/max(val,e); } );
    // Compute the Gradient
    Grad=c+A.t()*temp;
    // Compute the Hessian matrix
    vec dsq=temp;
    dsq.transform( [](double val) { return val*val; } );
    H=A.t()*diagmat(dsq)*A;

    vec dx = -solve(H, Grad, solve_opts::fast);
    double lambda = -as_scalar(Grad.t() * dx);
    if (lambda / 2 < threshold) {
      break;
    }
    //Backtracking Line Search
    f=ComputeF(A,b,c,p);
    double mu=1;
    while (ComputeF(A,b,c,p+mu*dx)>f+alphaparam*mu*as_scalar(Grad.t()*dx)){
      mu=mu*betaparam;
    }

    // update p
    p = p + mu*dx;
    iter++;
  }
}
void NewtonsMethoddual(mat &A, vec &b, vec &c, vec &p,double T,double& dualf,vec& dualArg,double& f) {
  int max_iter = 100;
  int iter = 0;
  int n=A.n_cols;
  int m=A.n_rows;
  vec Grad = vec(n);
  Grad.fill(1);
  mat H = mat(n, n);
  vec temp,dsq;
  while (iter < max_iter && norm(Grad) > threshold) {
    temp = b-A * p;
    temp.transform( [](double val) { return 1/max(val,e); } );
    // Compute the Gradient
    Grad=c+A.t()*temp;
    // Compute the Hessian matrix
    dsq=temp;
    dsq.transform( [](double val) { return val*val; } );
    H=A.t()*diagmat(dsq)*A;

    vec dx = -solve(H, Grad, solve_opts::fast);
    double lambda = -as_scalar(Grad.t() * dx);
    if (lambda / 2 < threshold) {
      break;
    }
    //Backtracking Line Search
    f=ComputeF(A,b,c,p);
    double mu=1;
    while (ComputeF(A,b,c,p+mu*dx)>f+alphaparam*mu*as_scalar(Grad.t()*dx)){
      mu=mu*betaparam;
    }
    // update p
    p = p + mu*dx;


    iter++;
  }
  f=ComputeF(A,b,c,p);
  dualf=f-m/T;
  dualArg=1/T*temp;

}
