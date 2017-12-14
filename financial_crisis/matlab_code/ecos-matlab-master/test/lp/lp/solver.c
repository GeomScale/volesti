#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[7] = {0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0};
  b[2] = p->b;
  static double c[8] = {1}; /* rest = {0.0}; */
  static double h[2]; /* = {0.0}; */
  static idxint q[0] = {};
  static idxint Gjc[9] = {0, 0, 1, 2, 2, 2, 2, 2, 2};
  static idxint Gir[2] = {0, 1};
  static double Gpr[2] = {-1.0, -1.0};
  static idxint Ajc[9] = {0, 1, 2, 3, 5, 7, 10, 13, 16};
  static idxint Air[16] = {6, 3, 4, 0, 1, 0, 2, 1, 3, 6, 1, 4, 6, 3, 4, 5};
  static double Apr[16] = {-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,0,1.0,0,0,1.0,0,-1.0,-1.0,1.0};
  static const idxint A_ind_map[2] = {7, 10};
  ptr = (double *) p->A;
  for (i = 0; i < 2; ++i) {
    Apr[ A_ind_map[i] ] = *ptr++;
  }
  static const idxint ct_ind_map[2] = {9, 12};
  ptr = (double *) p->ct;
  for (i = 0; i < 2; ++i) {
    Apr[ ct_ind_map[i] ] = *ptr++;
  }
  return ECOS_setup(8 /* num vars */, 2 /* num cone constraints */, 7 /* num eq constraints */, 2 /* num linear cones */, 0 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  for(i = 0; i < 2; ++i) {
    sol->x[i] = w->x[i + 5];
  }
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

