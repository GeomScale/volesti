#include "lp_lib.h"
#include "R.h"
#include "Rinternals.h"
#include "Rdefines.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Utils.h"

void R_init_lpSolveAPI(DllInfo *info);
lprec* lprecPointerFromSEXP(SEXP Slprec);
int __WINAPI RlpSolveAbortFunction(lprec *lp, void *userhandle);
void __WINAPI RlpSolveLogFunction(lprec *lp, void *userhandle, char *buf);
void RlpsHS(lprec *lp, unsigned char status);


