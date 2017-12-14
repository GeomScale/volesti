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
      p.A[i][j] = randu();
    }
  }
  p.b = randu();
  for(i = 0; i < 1; ++i) {
    for(j = 0; j < 2; ++j) {
      p.ct[i][j] = randu();
    }
  }
  pwork *w = setup(&p);
  int flag = 0;
  if(w!=NULL) {
    vars v;
    flag = solve(w, &v);
    cleanup(w);
  }
  return flag;
}
