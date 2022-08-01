#ifndef INTERIORPOINTMETHOD_H
#define INTERIORPOINTMETHOD_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

double IPM(mat& A,vec& b,vec& c,vec& p);
void NewtonsMethod(mat& A,vec& b,vec& c,vec& p);
void NewtonsMethoddual(mat &A, vec &b, vec &c, vec &p,double T,double& dualf,vec& dualArg,double& f);
double IPMprimaldual(mat &A, vec &b, vec &c, vec &p);

#endif