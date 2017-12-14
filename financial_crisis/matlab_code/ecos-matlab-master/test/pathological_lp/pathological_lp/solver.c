#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[6] = {0.0, 0.0, 0.0, 0.0, 0, 0};
  ptr = (double *) p->b;
  for(i = 4; i < 6; ++i) {
    b[i] = *ptr++;
  }
  static double c[7] = {1}; /* rest = {0.0}; */
  static double h[2]; /* = {0.0}; */
  static idxint q[0] = {};
  static idxint Gjc[8] = {0, 0, 1, 2, 2, 2, 2, 2};
  static idxint Gir[2] = {0, 1};
  static double Gpr[2] = {-1.0, -1.0};
  static idxint Ajc[8] = {0, 2, 3, 4, 6, 8, 10, 12};
  static idxint Air[12] = {2, 3, 0, 1, 0, 2, 1, 3, 0, 4, 1, 5};
  static double Apr[12] = {0,0,1.0,1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0,-1.0,1.0};
  static const idxint A_ind_map[2] = {0, 1};
  ptr = (double *) p->A;
  for (i = 0; i < 2; ++i) {
    Apr[ A_ind_map[i] ] = *ptr++;
  }
  return ECOS_setup(7 /* num vars */, 2 /* num cone constraints */, 6 /* num eq constraints */, 2 /* num linear cones */, 0 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  sol->x = w->x[0];
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

