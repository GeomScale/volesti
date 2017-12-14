#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0};
  b[5] = p->b;
  b[11] = p->b;
  static double c[20] = {-1}; /* rest = {0.0}; */
  static double h[8]; /* = {0.0}; */
  static idxint q[1] = {2};
  static idxint Gjc[21] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
  static idxint Gir[8] = {0, 6, 7, 1, 2, 3, 4, 5};
  static double Gpr[8] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  static idxint Ajc[21] = {0, 2, 3, 4, 7, 8, 9, 10, 11, 12, 14, 17, 20, 22, 26, 29, 31, 33, 35, 37, 39};
  static idxint Air[39] = {12, 13, 4, 4, 0, 3, 14, 6, 7, 8, 12, 13, 0, 1, 1, 9, 14, 2, 10, 16, 2, 3, 4, 5, 6, 11, 6, 7, 8, 7, 9, 8, 10, 12, 15, 13, 16, 14, 15};
  static double Apr[39] = {-1.0,-1.0,1.0,1.0,1.0,0,-1.0,1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0,1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0};
  static const idxint a_ind_map[2] = {5, 14};
  for (i = 0; i < 2; ++i) {
    Apr[ a_ind_map[i] ] = p->a;
  }
  return ECOS_setup(20 /* num vars */, 8 /* num cone constraints */, 17 /* num eq constraints */, 6 /* num linear cones */, 1 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  sol->x1 = w->x[10];
  sol->x2 = w->x[3];
  sol->x3 = w->x[11];
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

