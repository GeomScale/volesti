#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[17] = {0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5, -0.5, 0.0, 0.0, 0.0};
  b[2] = p->a;
  b[5] = p->b;
  static double c[23] = {1}; /* rest = {0.0}; */
  static double h[13]; /* = {0.0}; */
  static idxint q[2] = {3, 3};
  static idxint Gjc[24] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13};
  static idxint Gir[13] = {0, 1, 2, 3, 4, 5, 7, 8, 9, 6, 12, 10, 11};
  static double Gpr[13] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  static idxint Ajc[24] = {0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 21, 23, 26, 28, 31, 33, 36, 38, 40};
  static idxint Air[40] = {6, 0, 0, 1, 1, 3, 3, 4, 4, 9, 10, 11, 9, 10, 12, 13, 0, 2, 1, 8, 15, 3, 5, 4, 7, 16, 6, 7, 6, 9, 10, 7, 8, 12, 13, 14, 14, 15, 14, 16};
  static double Apr[40] = {-1.0,1.0,1.0,1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,0.5,0.5,-1.0,-1.0,-1.0,1.0,-1.0,-1.0,0,-1.0,1.0,1.0,1.0,0,1.0,-1.0,1.0,0.5,-0.5,1.0,-1.0,0.5,-0.5,-1.0,1.0,-1.0,-1.0,-1.0};
  Apr[20] = p->a;
  Apr[25] = p->b;
  return ECOS_setup(23 /* num vars */, 13 /* num cone constraints */, 17 /* num eq constraints */, 7 /* num linear cones */, 2 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  sol->x1 = w->x[14];
  sol->x2 = w->x[16];
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

