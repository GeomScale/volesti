// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#ifndef HEADER_lp_report
#define HEADER_lp_report

#ifdef __cplusplus
extern "C" {
#endif

/* General information functions */
char * __VACALL explain(lprec *lp, char *format, ...);
void __VACALL report(lprec *lp, int level, char *format, ...);

/* Prototypes for debugging and general data dumps */
void debug_print(lprec *lp, char *format, ...);
void debug_print_solution(lprec *lp);
void debug_print_bounds(lprec *lp, LPSREAL *upbo, LPSREAL *lowbo);
void blockWriteLREAL(FILE *output, char *label, LREAL *vector, int first, int last);
void blockWriteAMAT(FILE *output, const char *label, lprec* lp, int first, int last);
void blockWriteBMAT(FILE *output, const char *label, lprec* lp, int first, int last);


/* Model reporting headers */
void REPORT_objective(lprec *lp);
void REPORT_solution(lprec *lp, int columns);
void REPORT_constraints(lprec *lp, int columns);
void REPORT_duals(lprec *lp);
void REPORT_extended(lprec *lp);

/* Other rarely used, but sometimes extremely useful reports */
void REPORT_constraintinfo(lprec *lp, char *datainfo);
void REPORT_modelinfo(lprec *lp, MYBOOL doName, char *datainfo);
void REPORT_lp(lprec *lp);
MYBOOL REPORT_tableau(lprec *lp);
void REPORT_scales(lprec *lp);
MYBOOL REPORT_debugdump(lprec *lp, char *filename, MYBOOL livedata);
MYBOOL REPORT_mat_mmsave(lprec *lp, char *filename, int *colndx, MYBOOL includeOF, char *infotext);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_report */

