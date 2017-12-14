#include "solver.h"

pwork *setup(params *p)
{
  idxint i = 0;
  double *ptr;
  static double b[20] = {0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0};
  b[4] = p->b1;
  b[9] = p->b2;
  b[14] = p->b3;
  b[19] = p->b4;
  static double c[23] = {-1}; /* rest = {0.0}; */
  static double h[4]; /* = {0.0}; */
  static idxint q[0] = {};
  static idxint Gjc[24] = {0, 0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  static idxint Gir[4] = {0, 1, 2, 3};
  static double Gpr[4] = {-1.0, -1.0, -1.0, -1.0};
  static idxint Ajc[24] = {0, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48};
  static idxint Air[48] = {3, 8, 13, 18, 0, 5, 10, 15, 0, 1, 0, 4, 1, 2, 1, 3, 2, 7, 12, 17, 2, 7, 12, 17, 5, 6, 5, 9, 6, 7, 6, 8, 10, 11, 10, 14, 11, 12, 11, 13, 15, 16, 15, 19, 16, 17, 16, 18};
  static double Apr[48] = {0,0,0,0,1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,0,0,0,0,0,0,0,0,1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0};
  static const idxint a1t_ind_map[2] = {16, 20};
  ptr = (double *) p->a1t;
  for (i = 0; i < 2; ++i) {
    Apr[ a1t_ind_map[i] ] = *ptr++;
  }
  static const idxint a2t_ind_map[2] = {17, 21};
  ptr = (double *) p->a2t;
  for (i = 0; i < 2; ++i) {
    Apr[ a2t_ind_map[i] ] = *ptr++;
  }
  static const idxint a3t_ind_map[2] = {18, 22};
  ptr = (double *) p->a3t;
  for (i = 0; i < 2; ++i) {
    Apr[ a3t_ind_map[i] ] = *ptr++;
  }
  static const idxint a4t_ind_map[2] = {19, 23};
  ptr = (double *) p->a4t;
  for (i = 0; i < 2; ++i) {
    Apr[ a4t_ind_map[i] ] = *ptr++;
  }
  Apr[0] = p->na1;
  Apr[1] = p->na2;
  Apr[2] = p->na3;
  Apr[3] = p->na4;
  return ECOS_setup(23 /* num vars */, 4 /* num cone constraints */, 20 /* num eq constraints */, 4 /* num linear cones */, 0 /* num second-order cones */, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h ,b);
}

int solve(pwork *w, vars *sol)
{
  idxint i = 0;
  int exitflag = ECOS_solve(w);
  sol->r = w->x[0];
  for(i = 0; i < 2; ++i) {
    sol->x[i] = w->x[i + 9];
  }
  return exitflag;
}

void cleanup(pwork *w)
{
  ECOS_cleanup(w,0);
}

