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
  for(i = 0; i < 1; ++i) {
    for(j = 0; j < 2; ++j) {
      p.a1t[i][j] = randu();
    }
  }
  for(i = 0; i < 1; ++i) {
    for(j = 0; j < 2; ++j) {
      p.a2t[i][j] = randu();
    }
  }
  for(i = 0; i < 1; ++i) {
    for(j = 0; j < 2; ++j) {
      p.a3t[i][j] = randu();
    }
  }
  for(i = 0; i < 1; ++i) {
    for(j = 0; j < 2; ++j) {
      p.a4t[i][j] = randu();
    }
  }
  p.b1 = randu();
  p.b2 = randu();
  p.b3 = randu();
  p.b4 = randu();
  p.na1 = randu();
  p.na2 = randu();
  p.na3 = randu();
  p.na4 = randu();
  pwork *w = setup(&p);
  int flag = 0;
  if(w!=NULL) {
    vars v;
    flag = solve(w, &v);
    cleanup(w);
  }
  return flag;
}
