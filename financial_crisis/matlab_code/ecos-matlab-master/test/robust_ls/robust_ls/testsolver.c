#include "solver.h"
#include <stdlib.h> /* for random numbers */
#include <time.h> /* for seed */

double randu()
{
  return ((double) rand())/RAND_MAX;
}

int main(int argc, char **argv)
{
  srand( time(NULL) );
  idxint i = 0;
  idxint j = 0;
  params p;
  p.a11 = randu();
  p.a12 = randu();
  p.a21 = randu();
  p.a22 = randu();
  p.a31 = randu();
  p.a32 = randu();
  p.a41 = randu();
  p.a42 = randu();
  p.b1 = randu();
  p.b2 = randu();
  p.b3 = randu();
  p.b4 = randu();
  p.c11 = randu();
  p.c12 = randu();
  p.d = randu();
  p.e = randu();
  pwork *w = setup(&p);
  int flag = 0;
  if(w!=NULL) {
    vars v;
    flag = solve(w, &v);
    cleanup(w);
  }
  return flag;
}
