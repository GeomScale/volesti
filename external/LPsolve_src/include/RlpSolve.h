#include "lp_lib.h"
// #include "R.h"
// #include "Rinternals.h"
#define STRICT_R_HEADERS // To disable some R code else I get an compile error
#include <Rdefines.h>
// #include "R_ext/Rdynload.h"
// #include "R_ext/Utils.h"

void R_init_lpSolveAPI(DllInfo *info);
lprec* lprecPointerFromSEXP(SEXP Slprec);
int __WINAPI RlpSolveAbortFunction(lprec *lp, void *userhandle);
void __WINAPI RlpSolveLogFunction(lprec *lp, void *userhandle, char *buf);
void RlpsHS(lprec *lp, unsigned char status);


