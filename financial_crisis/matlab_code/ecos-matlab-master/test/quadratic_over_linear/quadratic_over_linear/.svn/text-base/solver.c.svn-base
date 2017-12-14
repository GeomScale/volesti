#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[14] = {0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5};
  b[1] = p->u;
  b[3] = p->l;
  b[5] = p->u;
  b[7] = p->l;
  static double c[16] = {1}; /* rest = {0.0}; */
  static double h[11]; /* = {0.0}; */
  static idxint q[2] = {3, 3};
  static idxint Gjc[17] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 11, 11, 11, 11, 11};
  static idxint Gir[11] = {0, 1, 2, 3, 5, 6, 7, 4, 10, 8, 9};
  static double Gpr[11] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  static idxint Ajc[17] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 16, 20, 24, 29, 32};
  static idxint Air[32] = {8, 0, 2, 4, 6, 9, 10, 11, 9, 10, 12, 13, 0, 2, 8, 11, 0, 1, 4, 5, 2, 3, 6, 7, 4, 6, 11, 12, 13, 8, 9, 10};
  static double Apr[32] = {-1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,0.5,0.5,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,0.5,-0.5,1.0,0.5,-0.5};

  return ECOS_setup(16 /* num vars */, 11 /* num cone constraints */, 14 /* num eq constraints */, 5 /* num linear cones */, 2 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  sol->x1 = w->x[11];
  sol->x2 = w->x[14];
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

