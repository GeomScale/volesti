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
  int i = 0;
  int j = 0;
  params p;
  p.B = 20; // randu();
  for(i = 0; i < 5; ++i) {
    for(j = 0; j < 5; ++j) {
      p.D[i][j] = 0.0; //randu();
    }
  }
  p.D[0][0] = 0.9115;
  p.D[1][1] = 0.7650;
  p.D[2][2] = 0.7414;
  p.D[3][3] = 0.9577;
  p.D[4][4] = 0.5346;

  for(i = 0; i < 5; ++i) {
    for(j = 0; j < 2; ++j) {
      p.F[i][j] = 0;
    }
  }
  p.F[0][0] = -0.6156; p.F[0][1] = -1.4023;
  p.F[1][0] = 0.7481; p.F[1][1] = -1.4224;
  p.F[2][0] = -0.1924; p.F[2][1] = 0.4882;
  p.F[3][0] = 0.8886; p.F[3][1] = -0.1774;
  p.F[4][0] = -0.7648; p.F[4][1] = -0.1961;

  p.gamma = 1;
  for(i = 0; i < 5; ++i) {
    p.mu[i] = 1;
  }
  pwork *w = setup(&p);
  int flag = 0;
  if(w!=NULL) {
    vars v;
    flag = solve(w, &v);
    for(i = 0; i < 5; i++){
        printf("%f\n", v.x[i]);
    }
    cleanup(w);
  }
  return flag;
}
