#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[8] = {0.0, 0, 0.0, 0, 0.0, 0.0, 0.0, 0.0};
  b[1] = p->a;
  b[3] = p->b;
  static double c[11] = {-1}; /* rest = {0.0}; */
  static double h[7]; /* = {0.0}; */
  static idxint q[1] = {3};
  static idxint Gjc[12] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7};
  static idxint Gir[7] = {0, 1, 2, 4, 5, 6, 3};
  static double Gpr[7] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  static idxint Ajc[12] = {0, 1, 2, 3, 4, 5, 6, 7, 10, 12, 18, 20};
  static idxint Air[20] = {5, 0, 2, 4, 6, 7, 5, 4, 6, 7, 0, 1, 0, 2, 4, 5, 6, 7, 2, 3};
  static double Apr[20] = {-1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,0.5,0.5,1.0,1.0,-1.0,1.0,-1.0,1.0,0.5,-0.5,-1.0,1.0};

  return ECOS_setup(11 /* num vars */, 7 /* num cone constraints */, 8 /* num eq constraints */, 4 /* num linear cones */, 1 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  sol->x1 = w->x[9];
  sol->x2 = w->x[7];
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

