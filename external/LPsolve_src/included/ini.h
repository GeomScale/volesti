// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#include <stdio.h>

#ifdef __cplusplus
__EXTERN_C {
#endif

extern FILE *ini_create(char *filename);
extern FILE *ini_open(char *filename);
extern void ini_writecomment(FILE *fp, char *comment);
extern void ini_writeheader(FILE *fp, char *header, int addnewline);
extern void ini_writedata(FILE *fp, char *name, char *data);
extern int ini_readdata(FILE *fp, char *data, int szdata, int withcomment);
extern void ini_close(FILE *fp);

#ifdef __cplusplus
}
#endif
