
/*
    Mixed integer programming optimization drivers for lp_solve v5.0+
   ----------------------------------------------------------------------------------
    Author:        Michel Berkelaar (to lp_solve v3.2),
                   Kjell Eikland
    Contact:
    License terms: LGPL.

    Requires:      stdarg.h, lp_lib.h

    Release notes:
    v5.0.0 3   1 January 2004      New unit isolating reporting routines.
    v5.2.0.0   1 December 2005     Addition of Matrix Market writing function.

   ----------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "lp_lib.h"
#include "lp_scale.h"
#include "commonlib.h"
#include "lp_report.h"

#include "mmio.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif

/* Define buffer-size controled function mapping */
# if defined _MSC_VER
#  define vsnprintf _vsnprintf
# endif

/* Various reporting functions for lp_solve                                  */
/* ------------------------------------------------------------------------- */

/* First define general utilties for reporting and output */
char * __VACALL explain(lprec *lp, char *format, ...)
{
  char buff[DEF_STRBUFSIZE+1];
  va_list ap;

  va_start(ap, format);
  vsnprintf(buff, DEF_STRBUFSIZE, format, ap);
  va_end(ap);
  allocCHAR(lp, &(lp->ex_status), (int) strlen(buff), AUTOMATIC);
  strcpy(lp->ex_status, buff);
  return( lp->ex_status );
}
void __VACALL report(lprec *lp, int level, char *format, ...)
{
  char buff[DEF_STRBUFSIZE+1];
  va_list ap;

  if(lp == NULL) {
    va_start(ap, format);
    //vfprintf(stderr, format, ap);
    va_end(ap);
  }
  else if(level <= lp->verbose) {
    if(lp->writelog != NULL) {
      va_start(ap, format);
      vsnprintf(buff, DEF_STRBUFSIZE, format, ap);
      va_end(ap);
      lp->writelog(lp, lp->loghandle, buff);
    }
    if(lp->outstream != NULL) {
      va_start(ap, format);
      //vfprintf(lp->outstream, format, ap);
      va_end(ap);
      //if(lp->outstream != stdout)
      //  fflush(lp->outstream);
    }
  }
#ifdef xParanoia
  if(level == CRITICAL)
    raise(SIGSEGV);
#endif
}

STATIC void print_indent(lprec *lp)
{
  int i;

  report(lp, NEUTRAL, "%2d", lp->bb_level);
  if(lp->bb_level < 50) /* useless otherwise */
    for(i = lp->bb_level; i > 0; i--)
      report(lp, NEUTRAL, "--");
  else
    report(lp, NEUTRAL, " *** too deep ***");
  report(lp, NEUTRAL, "> ");
} /* print_indent */

STATIC void debug_print(lprec *lp, char *format, ...)
{
  va_list ap;

  if(lp->bb_trace) {
    print_indent(lp);
    if (lp == NULL)
    {
      va_start(ap, format);
      //vfprintf(stderr, format, ap);
      va_end(ap);
      //fputc('\n', stderr);
    }
    else if(lp->debuginfo != NULL)
    {
      char buff[DEF_STRBUFSIZE+1];
      va_start(ap, format);
      vsnprintf(buff, DEF_STRBUFSIZE, format, ap);
      va_end(ap);
      lp->debuginfo(lp, lp->loghandle, buff);
    }
  }
} /* debug_print */

STATIC void debug_print_solution(lprec *lp)
{
  int i;

  if(lp->bb_trace)
    for (i = lp->rows + 1; i <= lp->sum; i++) {
      print_indent(lp);
      report(lp, NEUTRAL, "%s " RESULTVALUEMASK "\n",
                 get_col_name(lp, i - lp->rows),
                (double)lp->solution[i]);
    }
} /* debug_print_solution */

STATIC void debug_print_bounds(lprec *lp, LPSREAL *upbo, LPSREAL *lowbo)
{
  int i;

  if(lp->bb_trace)
    for(i = lp->rows + 1; i <= lp->sum; i++) {
      if(lowbo[i] == upbo[i]) {
        print_indent(lp);
        report(lp, NEUTRAL, "%s = " RESULTVALUEMASK "\n", get_col_name(lp, i - lp->rows),
                             (double)lowbo[i]);
      }
      else {
        if(lowbo[i] != 0) {
          print_indent(lp);
          report(lp, NEUTRAL, "%s > " RESULTVALUEMASK "\n", get_col_name(lp, i - lp->rows),
                               (double)lowbo[i]);
        }
        if(upbo[i] != lp->infinite) {
          print_indent(lp);
          report(lp, NEUTRAL, "%s < " RESULTVALUEMASK "\n", get_col_name(lp, i - lp->rows),
                               (double)upbo[i]);
    }
      }
    }
} /* debug_print_bounds */

/* List a vector of LREAL values for the given index range */
void blockWriteLREAL(FILE *output, char *label, LREAL *vector, int first, int last)
{
  int i, k = 0;

  //fprintf(output, "%s", label);
  //fprintf(output, "\n");
  for(i = first; i <= last; i++) {
    //fprintf(output, " %18g", vector[i]);
    k++;
    if(my_mod(k, 4) == 0) {
      //fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 4) != 0)
    ;//fprintf(output, "\n");
}

/* List the current user data matrix columns over the selected row range */
void blockWriteAMAT(FILE *output, const char *label, lprec* lp, int first, int last)
{
  int    i, j, k = 0;
  int    nzb, nze, jb;
  double hold;
  MATrec *mat = lp->matA;

  if(!mat_validate(mat))
    return;
  if(first < 0)
    first = 0;
  if(last < 0)
    last = lp->rows;

  //fprintf(output, "%s", label);
  //fprintf(output, "\n");

  if(first == 0) {
    for(j = 1; j <= lp->columns; j++) {
      hold = get_mat(lp, 0, j);
      //fprintf(output, " %18g", hold);
      k++;
      if(my_mod(k, 4) == 0) {
        //fprintf(output, "\n");
        k = 0;
      }
    }
    if(my_mod(k, 4) != 0) {
      //fprintf(output, "\n");
      k = 0;
    }
    first++;
  }
  nze = mat->row_end[first-1];
  for(i = first; i <= last; i++) {
    nzb = nze;
    nze = mat->row_end[i];
    if(nzb >= nze)
      jb = lp->columns+1;
    else
      jb = ROW_MAT_COLNR(nzb);
    for(j = 1; j <= lp->columns; j++) {
      if(j < jb)
        hold = 0;
      else {
        hold = get_mat(lp, i, j);
        nzb++;
        if(nzb < nze)
          jb = ROW_MAT_COLNR(nzb);
        else
          jb = lp->columns+1;
      }
      //fprintf(output, " %18g", hold);
      k++;
      if(my_mod(k, 4) == 0) {
        //fprintf(output, "\n");
        k = 0;
      }
    }
    if(my_mod(k, 4) != 0) {
      //fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 4) != 0)
    ;//fprintf(output, "\n");
}

/* List the current basis matrix columns over the selected row range */
void blockWriteBMAT(FILE *output, const char *label, lprec* lp, int first, int last)
{
  int    i, j, jb, k = 0;
  double hold;

  if(first < 0)
    first = 0;
  if(last < 0)
    last = lp->rows;

  //fprintf(output, "%s", label);
  //fprintf(output, "\n");

  for(i = first; i <= last; i++) {
    for(j = 1; j <= lp->rows; j++) {
      jb = lp->var_basic[j];
      if(jb <= lp->rows) {
        if(jb == i)
          hold = 1;
        else
          hold = 0;
      }
      else
        hold = get_mat(lp, i, j);
      if(i == 0)
        modifyOF1(lp, jb, &hold, 1);
      hold = unscaled_mat(lp, hold, i, jb);
      ;//fprintf(output, " %18g", hold);
      k++;
      if(my_mod(k, 4) == 0) {
        ;//fprintf(output, "\n");
        k = 0;
      }
    }
    if(my_mod(k, 4) != 0) {
      ;//fprintf(output, "\n");
      k = 0;
    }
  }
  if(my_mod(k, 4) != 0)
    ;//fprintf(output, "\n");
}

/* Do a generic readable data dump of key lp_solve model variables;
   principally for run difference and debugging purposes */
MYBOOL REPORT_debugdump(lprec *lp, char *filename, MYBOOL livedata)
{
  return 1;
}


/* High level reports for model results */

void REPORT_objective(lprec *lp)
{}

void REPORT_solution(lprec *lp, int columns)
{} /* REPORT_solution */

void REPORT_constraints(lprec *lp, int columns)
{}

void REPORT_duals(lprec *lp)
{
  int i;
  LPSREAL *duals, *dualsfrom, *dualstill, *objfrom, *objtill, *objfromvalue;
  MYBOOL ret;

  if(lp->outstream == NULL)
    return;

  ret = get_ptr_sensitivity_objex(lp, &objfrom, &objtill, &objfromvalue, NULL);
  /*
  if(ret) {
    //fprintf(lp->outstream, "\nObjective function limits:\n");
    //fprintf(lp->outstream, "                                 From            Till       FromValue\n");
    //for(i = 1; i <= lp->columns; i++)
    //  if(!is_splitvar(lp, i))
        //fprintf(lp->outstream, "%-20s  %15.7g %15.7g %15.7g\n", get_col_name(lp, i),
        // (double)objfrom[i - 1], (double)objtill[i - 1], (double)objfromvalue[i - 1]);
  }*/

  ret = get_ptr_sensitivity_rhs(lp, &duals, &dualsfrom, &dualstill);
  /*
  if(ret) {
    fprintf(lp->outstream, "\nDual values with from - till limits:\n");
    fprintf(lp->outstream, "                           Dual value            From            Till\n");
    for(i = 1; i <= lp->sum; i++)
      fprintf(lp->outstream, "%-20s  %15.7g %15.7g %15.7g\n",
              (i <= lp->rows) ? get_row_name(lp, i) : get_col_name(lp, i - lp->rows),
              (double)duals[i - 1], (double)dualsfrom[i - 1], (double)dualstill[i - 1]);
    fflush(lp->outstream);
  }*/
}

/* Printing of sensitivity analysis reports */
void REPORT_extended(lprec *lp)
{
  int  i, j;
  LPSREAL hold;
  LPSREAL *duals, *dualsfrom, *dualstill, *objfrom, *objtill;
  MYBOOL ret;

  ret = get_ptr_sensitivity_obj(lp, &objfrom, &objtill);
  report(lp, NORMAL, " \n");
  report(lp, NORMAL, "Primal objective:\n");
  report(lp, NORMAL, " \n");
  report(lp, NORMAL, "  Column name                      Value   Objective         Min         Max\n");
  report(lp, NORMAL, "  --------------------------------------------------------------------------\n");
  for(j = 1; j <= lp->columns; j++) {
    hold = get_mat(lp,0,j);
    report(lp, NORMAL, "  %-25s " MPSVALUEMASK MPSVALUEMASK MPSVALUEMASK MPSVALUEMASK "\n",
           get_col_name(lp,j),
           my_precision(hold,lp->epsprimal),
           my_precision(hold*lp->best_solution[lp->rows+j],lp->epsprimal),
           my_precision((ret) ? objfrom[j - 1] : 0.0,lp->epsprimal),
           my_precision((ret) ? objtill[j - 1] : 0.0,lp->epsprimal));
  }
  report(lp, NORMAL, " \n");

  ret = get_ptr_sensitivity_rhs(lp, &duals, &dualsfrom, &dualstill);
  report(lp, NORMAL, "Primal variables:\n");
  report(lp, NORMAL, " \n");
  report(lp, NORMAL, "  Column name                      Value       Slack         Min         Max\n");
  report(lp, NORMAL, "  --------------------------------------------------------------------------\n");
  for(j = 1; j <= lp->columns; j++)
    report(lp, NORMAL, "  %-25s " MPSVALUEMASK MPSVALUEMASK MPSVALUEMASK MPSVALUEMASK "\n",
           get_col_name(lp,j),
           my_precision(lp->best_solution[lp->rows+j],lp->epsprimal),
           my_precision(my_inflimit(lp, (ret) ? duals[lp->rows+j-1] : 0.0),lp->epsprimal),
           my_precision((ret) ? dualsfrom[lp->rows+j-1] : 0.0,lp->epsprimal),
           my_precision((ret) ? dualstill[lp->rows+j-1] : 0.0,lp->epsprimal));

  report(lp, NORMAL, " \n");
  report(lp, NORMAL, "Dual variables:\n");
  report(lp, NORMAL, " \n");
  report(lp, NORMAL, "  Row name                         Value       Slack         Min         Max\n");
  report(lp, NORMAL, "  --------------------------------------------------------------------------\n");
  for(i = 1; i <= lp->rows; i++)
    report(lp, NORMAL, "  %-25s " MPSVALUEMASK MPSVALUEMASK MPSVALUEMASK MPSVALUEMASK "\n",
           get_row_name(lp,i),
           my_precision((ret) ? duals[i - 1] : 0.0, lp->epsprimal),
           my_precision(lp->best_solution[i], lp->epsprimal),
           my_precision((ret) ? dualsfrom[i - 1] : 0.0,lp->epsprimal),
           my_precision((ret) ? dualstill[i - 1] : 0.0,lp->epsprimal));

  report(lp, NORMAL, " \n");
}

/* A more readable lp-format report of the model; antiquated and not updated */
void REPORT_lp(lprec *lp)
{}

/* Report the scaling factors used; extremely rarely used */
void REPORT_scales(lprec *lp)
{
  int i, colMax;

  colMax = lp->columns;

  if(lp->outstream == NULL)
    return;
  /*
  if(lp->scaling_used) {
    fprintf(lp->outstream, "\nScale factors:\n");
    for(i = 0; i <= lp->rows + colMax; i++)
      fprintf(lp->outstream, "%-20s scaled at %g\n",
              (i <= lp->rows) ? get_row_name(lp, i) : get_col_name(lp, i - lp->rows),
        (double)lp->scalars[i]);
  }
  fflush(lp->outstream);
  */
}

/* Report the traditional tableau corresponding to the current basis */
MYBOOL REPORT_tableau(lprec *lp)
{
  int  j, row_nr, *coltarget;
  LPSREAL *prow = NULL;
  FILE *stream;// = lp->outstream;

  if(lp->outstream == NULL)
    return(FALSE);

  if(!lp->model_is_valid || !has_BFP(lp) ||
     (get_total_iter(lp) == 0) || (lp->spx_status == NOTRUN)) {
    lp->spx_status = NOTRUN;
    return(FALSE);
  }
  if(!allocREAL(lp, &prow,lp->sum + 1, TRUE)) {
    lp->spx_status = NOMEMORY;
    return(FALSE);
  }

  coltarget = (int *) mempool_obtainVector(lp->workarrays, lp->columns+1, sizeof(*coltarget));
  if(!get_colIndexA(lp, SCAN_USERVARS+USE_NONBASICVARS, coltarget, FALSE)) {
    mempool_releaseVector(lp->workarrays, (char *) coltarget, FALSE);
    return(FALSE);
  }
  for(row_nr = 1; (row_nr <= lp->rows + 1); row_nr++) {
    bsolve(lp, row_nr <= lp->rows ? row_nr : 0, prow, NULL, lp->epsmachine*DOUBLEROUND, 1.0);
    prod_xA(lp, coltarget, prow, NULL, lp->epsmachine, 1.0,
                                       prow, NULL, MAT_ROUNDDEFAULT);

  }

  mempool_releaseVector(lp->workarrays, (char *) coltarget, FALSE);
  FREE(prow);
  return(TRUE);
}

void REPORT_constraintinfo(lprec *lp, char *datainfo)
{
  int i, tally[ROWCLASS_MAX+1];

  MEMCLEAR(tally, ROWCLASS_MAX+1);
  for(i = 1; i <= lp->rows; i++)
    tally[get_constr_class(lp, i)]++;

  if(datainfo != NULL)
    report(lp, NORMAL, "%s\n", datainfo);

  for(i = 0; i <= ROWCLASS_MAX; i++)
    if(tally[i] > 0)
      report(lp, NORMAL, "%-15s %4d\n", get_str_constr_class(lp, i), tally[i]);
}

void REPORT_modelinfo(lprec *lp, MYBOOL doName, char *datainfo)
{
  if(doName) {
    report(lp, NORMAL, "\nModel name:  '%s' - run #%-5d\n",
                       get_lp_name(lp), lp->solvecount);
    report(lp, NORMAL, "Objective:   %simize(%s)\n",
                       my_if(is_maxim(lp), "Max", "Min"), get_row_name(lp, 0));
    report(lp, NORMAL, " \n");
  }
  if(datainfo != NULL)
    report(lp, NORMAL, "%s\n", datainfo);

  report(lp, NORMAL, "Model size:  %7d constraints, %7d variables, %12d non-zeros.\n",
         lp->rows, lp->columns, get_nonzeros(lp));
  if(GUB_count(lp)+SOS_count(lp) > 0)
  report(lp, NORMAL, "Var-types:   %7d integer,     %7d semi-cont.,     %7d SOS.\n",
         lp->int_vars, lp->sc_vars, lp->sos_vars);
  report(lp, NORMAL, "Sets:                             %7d GUB,            %7d SOS.\n",
                         GUB_count(lp), SOS_count(lp));
}

/* Save a matrix column subset to a MatrixMarket formatted file,
   say to export the basis matrix for further numerical analysis.
   If colndx is NULL, then the full constraint matrix is assumed. */
MYBOOL REPORT_mat_mmsave(lprec *lp, char *filename, int *colndx, MYBOOL includeOF, char *infotext)
{
return 1;
}

