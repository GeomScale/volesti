#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[19] = {0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0.0, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0};
  b[3] = p->l;
  b[7] = p->u;
  ptr = (double *) p->b;
  for(i = 15; i < 17; ++i) {
    b[i] = *ptr++;
  }
  static double c[27] = {1}; /* rest = {0.0}; */
  static double h[18]; /* = {0.0}; */
  static idxint q[5] = {2, 2, 2, 3, 3};
  static idxint Gjc[28] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18};
  static idxint Gir[18] = {0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 6, 8, 10, 7, 9, 11};
  static double Gpr[18] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  static idxint Ajc[28] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 10, 11, 12, 13, 14, 18, 22, 26, 30, 34, 37, 39, 41, 43, 45, 47, 49};
  static idxint Air[49] = {8, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12, 18, 18, 18, 0, 4, 13, 14, 1, 5, 13, 14, 2, 6, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 17, 11, 13, 12, 14, 11, 15, 12, 16, 17, 18};
  static double Apr[49] = {-1.0,1.0,1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,1.0,0,0,-1.0,1.0,0,0,-1.0,1.0,0,0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,0.5,-0.5,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0,-1.0,1.0,0,-1.0};
  static const idxint A_ind_map[6] = {16, 20, 24, 17, 21, 25};
  ptr = (double *) p->A;
  for (i = 0; i < 6; ++i) {
    Apr[ A_ind_map[i] ] = *ptr++;
  }
  Apr[47] = p->lambda;
  return ECOS_setup(27 /* num vars */, 18 /* num cone constraints */, 19 /* num eq constraints */, 6 /* num linear cones */, 5 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  for(i = 0; i < 3; ++i) {
    sol->x[i] = w->x[i + 15];
  }
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

