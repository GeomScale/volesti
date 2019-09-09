// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#include "lp_lib.h"
#include "lp_report.h"

#define MALLOC(ptr, nr, type)\
  ((((nr) == 0) || ((ptr = (type *) malloc((size_t)((nr) * sizeof(*ptr)))) == NULL)) ? \
   report(NULL, CRITICAL, "malloc of %d bytes failed on line %d of file %s\n",\
           (nr) * sizeof(*ptr), __LINE__, __FILE__), (ptr = NULL /* (void *) 0 */) : \
   ptr\
  )

#define CALLOC(ptr, nr, type)\
  ((((nr) == 0) || ((ptr = (type *) calloc((size_t)(nr), sizeof(*ptr))) == NULL)) ? \
   report(NULL, CRITICAL, "calloc of %d bytes failed on line %d of file %s\n",\
           (nr) * sizeof(*ptr), __LINE__, __FILE__), (ptr = NULL /* (void *) 0 */) : \
   ptr\
  )

#define REALLOC(ptr, nr, type)\
  ((((nr) == 0) || ((ptr = (type *) realloc(ptr, (size_t)((nr) * sizeof(*ptr)))) == NULL)) ? \
   report(NULL, CRITICAL, "realloc of %d bytes failed on line %d of file %s\n",\
           (nr) * sizeof(*ptr), __LINE__, __FILE__), (ptr = NULL /* (void *) 0 */) : \
   ptr\
  )

#if defined FREE
# undef FREE
#endif

#define FREE(ptr) if (ptr != NULL) {free(ptr), ptr = NULL;} else

#define MALLOCCPY(nptr, optr, nr, type)\
  (MALLOC(nptr, nr, type), (nptr != NULL) ? memcpy(nptr, optr, (size_t)((nr) * sizeof(*optr))) : 0)
