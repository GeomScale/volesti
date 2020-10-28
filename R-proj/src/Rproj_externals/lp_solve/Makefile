LP_SOLVE_CPPFLAGS=$(CPPFLAGS) -I../../../../external/LPsolve_src/include \
                  -I$(R_INCLUDE_DIR) \
                  -DYY_NEVER_INTERACTIVE \
                  -DPARSER_LP -DINVERSE_ACTIVE=INVERSE_LUSOL \
                  -DRoleIsExternalInvEngine

LP_SOLVE_SOURCES=colamd.c lp_MDO.c lp_mipbb.c lp_rlp.c mmio.c commonlib.c \
                 lp_MPS.c lp_params.c lp_scale.c lp_SOS.c \
                 lp_presolve.c lp_simplex.c yacc_read.c ini.c lp_crash.c \
                 lp_price.c lp_utils.c lp_Hash.c lp_lib.c lp_pricePSE.c \
                 lp_wlp.c lp_LUSOL.c lp_matrix.c lp_report.c lusol.c \
								 myblas.c

LP_SOLVE_OBJECTS=$(LP_SOLVE_SOURCES:.c=.o)

liblp_solve.a: $(LP_SOLVE_OBJECTS)
	$(AR) rc liblp_solve.a $(LP_SOLVE_OBJECTS) && $(RANLIB) liblp_solve.a

.c.o:
	$(CC) $(CFLAGS) $(CPICFLAGS) $(LP_SOLVE_CPPFLAGS) -c $< -o $@

clean:
	rm -rf $(LP_SOLVE_OBJECTS) liblp_solve.a


