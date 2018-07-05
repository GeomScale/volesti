#ifndef HEADER_lp_utils
#define HEADER_lp_utils

#ifdef FORTIFY

#include "lp_fortify.h"

#define allocCHAR allocCHAR_FORTIFY
#define allocMYBOOL allocMYBOOL_FORTIFY
#define allocINT allocINT_FORTIFY
#define allocREAL allocREAL_FORTIFY
#define allocLREAL allocLREAL_FORTIFY

#endif

#include "lp_types.h"

/* Temporary data storage arrays */
typedef struct _workarraysrec
{
  lprec     *lp;
  int       size;
  int       count;
  char      **vectorarray;
  int       *vectorsize;
} workarraysrec;

typedef struct _LLrec
{
  int       size;               /* The allocated list size */
  int       count;              /* The current entry count */
  int       firstitem;
  int       lastitem;
  int       *map;               /* The list of forward and backward-mapped entries */
} LLrec;

typedef struct _PVrec
{
  int       count;              /* The allocated list item count */
  int       *startpos;          /* Starting index of the current value */
  LPSREAL      *value;             /* The list of forward and backward-mapped entries */
  struct   _PVrec *parent;     /* The parent record in a pushed chain */
} PVrec;


#ifdef __cplusplus
extern "C" {
#endif

/* Put function headers here */
STATIC MYBOOL allocCHAR(lprec *lp, char **ptr, int size, MYBOOL clear);
STATIC MYBOOL allocMYBOOL(lprec *lp, MYBOOL **ptr, int size, MYBOOL clear);
STATIC MYBOOL allocINT(lprec *lp, int **ptr, int size, MYBOOL clear);
STATIC MYBOOL allocREAL(lprec *lp, LPSREAL **ptr, int size, MYBOOL clear);
STATIC MYBOOL allocLREAL(lprec *lp, LREAL **ptr, int size, MYBOOL clear);
STATIC MYBOOL allocFREE(lprec *lp, void **ptr);
LPSREAL *cloneREAL(lprec *lp, LPSREAL *origlist, int size);
MYBOOL *cloneMYBOOL(lprec *lp, MYBOOL *origlist, int size);
int *cloneINT(lprec *lp, int *origlist, int size);

int comp_bits(MYBOOL *bitarray1, MYBOOL *bitarray2, int items);

STATIC workarraysrec *mempool_create(lprec *lp);
STATIC char *mempool_obtainVector(workarraysrec *mempool, int count, int unitsize);
STATIC MYBOOL mempool_releaseVector(workarraysrec *mempool, char *memvector, MYBOOL forcefree);
STATIC MYBOOL mempool_free(workarraysrec **mempool);

STATIC void roundVector(LREAL *myvector, int endpos, LREAL roundzero);
STATIC LPSREAL normalizeVector(LPSREAL *myvector, int endpos);

STATIC void swapINT(int *item1, int *item2);
STATIC void swapREAL(LPSREAL *item1, LPSREAL *item2);
STATIC void swapPTR(void **item1, void **item2);
STATIC LPSREAL restoreINT(LPSREAL valREAL, LPSREAL epsilon);
STATIC LPSREAL roundToPrecision(LPSREAL value, LPSREAL precision);

STATIC int searchFor(int target, int *attributes, int size, int offset, MYBOOL absolute);

STATIC MYBOOL isINT(lprec *lp, LPSREAL value);
STATIC MYBOOL isOrigFixed(lprec *lp, int varno);
STATIC void chsign_bounds(LPSREAL *lobound, LPSREAL *upbound);
STATIC LPSREAL rand_uniform(lprec *lp, LPSREAL range);

/* Doubly linked list routines */
STATIC int createLink(int size, LLrec **linkmap, MYBOOL *usedpos);
STATIC MYBOOL freeLink(LLrec **linkmap);
STATIC int sizeLink(LLrec *linkmap);
STATIC MYBOOL isActiveLink(LLrec *linkmap, int itemnr);
STATIC int countActiveLink(LLrec *linkmap);
STATIC int countInactiveLink(LLrec *linkmap);
STATIC int firstActiveLink(LLrec *linkmap);
STATIC int lastActiveLink(LLrec *linkmap);
STATIC MYBOOL appendLink(LLrec *linkmap, int newitem);
STATIC MYBOOL insertLink(LLrec *linkmap, int afteritem, int newitem);
STATIC MYBOOL setLink(LLrec *linkmap, int newitem);
STATIC MYBOOL fillLink(LLrec *linkmap);
STATIC int nextActiveLink(LLrec *linkmap, int backitemnr);
STATIC int prevActiveLink(LLrec *linkmap, int forwitemnr);
STATIC int firstInactiveLink(LLrec *linkmap);
STATIC int lastInactiveLink(LLrec *linkmap);
STATIC int nextInactiveLink(LLrec *linkmap, int backitemnr);
STATIC int prevInactiveLink(LLrec *linkmap, int forwitemnr);
STATIC int removeLink(LLrec *linkmap, int itemnr);
STATIC LLrec *cloneLink(LLrec *sourcemap, int newsize, MYBOOL freesource);
STATIC int compareLink(LLrec *linkmap1, LLrec *linkmap2);
STATIC MYBOOL verifyLink(LLrec *linkmap, int itemnr, MYBOOL doappend);

/* Packed vector routines */
STATIC PVrec  *createPackedVector(int size, LPSREAL *values, int *workvector);
STATIC void   pushPackedVector(PVrec *PV, PVrec *parent);
STATIC MYBOOL unpackPackedVector(PVrec *PV, LPSREAL **target);
STATIC LPSREAL   getvaluePackedVector(PVrec *PV, int index);
STATIC PVrec  *popPackedVector(PVrec *PV);
STATIC MYBOOL freePackedVector(PVrec **PV);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_utils */

#ifdef FORTIFY

#if defined CODE_lp_utils && !defined CODE_lp_utils_
int _Fortify_ret;
#else
extern int _Fortify_ret;
#endif

#ifdef CODE_lp_utils
#define CODE_lp_utils_
#else
# undef allocCHAR
# undef allocMYBOOL
# undef allocINT
# undef allocREAL
# undef allocLREAL
# define allocCHAR(lp, ptr, size, clear) (Fortify_LINE(__LINE__), Fortify_FILE(__FILE__), _Fortify_ret = allocCHAR_FORTIFY(lp, ptr, size, clear), Fortify_LINE(0), Fortify_FILE(NULL), _Fortify_ret)
# define allocMYBOOL(lp, ptr, size, clear) (Fortify_LINE(__LINE__), Fortify_FILE(__FILE__), _Fortify_ret = allocMYBOOL_FORTIFY(lp, ptr, size, clear), Fortify_LINE(0), Fortify_FILE(NULL), _Fortify_ret)
# define allocINT(lp, ptr, size, clear) (Fortify_LINE(__LINE__), Fortify_FILE(__FILE__), _Fortify_ret = allocINT_FORTIFY(lp, ptr, size, clear), Fortify_LINE(0), Fortify_FILE(NULL), _Fortify_ret)
# define allocREAL(lp, ptr, size, clear) (Fortify_LINE(__LINE__), Fortify_FILE(__FILE__), _Fortify_ret = allocREAL_FORTIFY(lp, ptr, size, clear), Fortify_LINE(0), Fortify_FILE(NULL), _Fortify_ret)
# define allocLREAL(lp, ptr, size, clear) (Fortify_LINE(__LINE__), Fortify_FILE(__FILE__), _Fortify_ret = allocLREAL_FORTIFY(lp, ptr, size, clear), Fortify_LINE(0), Fortify_FILE(NULL), _Fortify_ret)
#endif

#endif

