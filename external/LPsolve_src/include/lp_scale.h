#ifndef HEADER_lp_scale
#define HEADER_lp_scale

#include "lp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Put function headers here */
STATIC MYBOOL scale_updatecolumns(lprec *lp, LPSREAL *scalechange, MYBOOL updateonly);
STATIC MYBOOL scale_updaterows(lprec *lp, LPSREAL *scalechange, MYBOOL updateonly);
STATIC MYBOOL scale_rows(lprec *lp, LPSREAL *scaledelta);
STATIC MYBOOL scale_columns(lprec *lp, LPSREAL *scaledelta);
STATIC void unscale_columns(lprec *lp);
STATIC LPSREAL scale(lprec *lp, LPSREAL *scaledelta);
STATIC LPSREAL scaled_mat(lprec *lp, LPSREAL value, int rownr, int colnr);
STATIC LPSREAL unscaled_mat(lprec *lp, LPSREAL value, int rownr, int colnr);
STATIC LPSREAL scaled_value(lprec *lp, LPSREAL value, int index);
STATIC LPSREAL unscaled_value(lprec *lp, LPSREAL value, int index);
STATIC MYBOOL scaleCR(lprec *lp, LPSREAL *scaledelta);
STATIC MYBOOL finalize_scaling(lprec *lp, LPSREAL *scaledelta);
STATIC LPSREAL auto_scale(lprec *lp);
void undoscale(lprec *lp);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_scale */

