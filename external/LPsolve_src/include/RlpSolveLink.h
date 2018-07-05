/************************************************************************************
  * The lp_solve API 
************************************************************************************/

/*******************************
  * Create/destroy model
*******************************/

SEXP RlpSolve_make_lp(SEXP Srows, SEXP Scolumns);
SEXP RlpSolve_copy_lp(SEXP Slp);
/*read_lp*/
SEXP RlpSolve_read_LP(SEXP Sfilename, SEXP Sverbose);
/*read_mps*/
/*read_freemps*/
SEXP RlpSolve_read_MPS(SEXP Sfilename, SEXP Soptions);
SEXP RlpSolve_read_freeMPS(SEXP Sfilename, SEXP Soptions);
/*SEXP RlpSolve_read_XLI(SEXP Sxliname, SEXP Smodelname, SEXP Sdataname,
                       SEXP Soptions);*/
SEXP RlpSolve_delete_lp(SEXP Slp);
/*free_lp*/

/*******************************
  * Build model
*******************************/

/*add_column*/
SEXP RlpSolve_add_columnex(SEXP Slp, SEXP Scolumn, SEXP Srowno);
/*str_add_column*/
SEXP RlpSolve_set_columnex(SEXP Slp, SEXP Scol_no, SEXP Scolumn, SEXP Srowno);
/*get_column*/
SEXP RlpSolve_get_columnex(SEXP Slp, SEXP Scol_nr);
/*add_constraint*/
SEXP RlpSolve_add_constraintex(SEXP Slp, SEXP Srow, SEXP Scolno,
                               SEXP Sconstr_type, SEXP Srh);
/*str_add_constraint*/
/*set_row*/
SEXP RlpSolve_set_rowex(SEXP Slp, SEXP Srow_no, SEXP Srow, SEXP Scolno);
SEXP RlpSolve_add_lag_con(SEXP Slp, SEXP Srow, SEXP Sconstr_type, SEXP Srh);
/*str_add_lag_con*/
SEXP RlpSolve_add_SOS(SEXP Slp, SEXP Sname, SEXP Ssostype, SEXP Spriority,
                      SEXP Ssosvars, SEXP Sweights);
SEXP RlpSolve_is_SOS_var(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_del_column(SEXP Slp, SEXP Scolumn);
SEXP RlpSolve_del_constraint(SEXP Slp, SEXP Sdel_row);
/*get_row*/
SEXP RlpSolve_get_rowex(SEXP Slp, SEXP Srow_nr);
SEXP RlpSolve_get_nameindex(SEXP Slp, SEXP Snames, SEXP Sisrow);
SEXP RlpSolve_is_infinite(SEXP Slp, SEXP Svalues);
SEXP RlpSolve_is_negative(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_resize_lp(SEXP Slp, SEXP Srows, SEXP Scolumns);
SEXP RlpSolve_set_add_rowmode(SEXP Slp, SEXP Sturnon);
SEXP RlpSolve_is_add_rowmode(SEXP Slp);
SEXP RlpSolve_set_binary(SEXP Slp, SEXP Scolumns, SEXP Smust_be_bin);
SEXP RlpSolve_is_binary(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_bounds(SEXP Slp, SEXP Scolumns, SEXP Slower, SEXP Supper);
SEXP RlpSolve_set_bounds_tighter(SEXP Slp, SEXP Stighten);
SEXP RlpSolve_get_bounds_tighter(SEXP Slp);
SEXP RlpSolve_set_col_names(SEXP Slp, SEXP Scolumns, SEXP Snames);
SEXP RlpSolve_get_col_names(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_get_origcol_names(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_constr_type(SEXP Slp, SEXP Srows, SEXP Scon_types);
SEXP RlpSolve_get_constr_type(SEXP Slp, SEXP Srows);
SEXP RlpSolve_is_constr_type(SEXP Slp, SEXP Srows, SEXP Smasks);
SEXP RlpSolve_set_unbounded(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_is_unbounded(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_infinite(SEXP Slp, SEXP Sinfinite);
SEXP RlpSolve_get_infinite(SEXP Slp);
SEXP RlpSolve_set_int(SEXP Slp, SEXP Scolumns, SEXP Smust_be_int);
SEXP RlpSolve_is_int(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_lowbo(SEXP Slp, SEXP Scolumns, SEXP Svalues);
SEXP RlpSolve_get_lowbo(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_lp_name(SEXP Slp, SEXP Slpname);
SEXP RlpSolve_get_lp_name(SEXP Slp);
SEXP RlpSolve_set_mat(SEXP Slp, SEXP Srow, SEXP Scolumn, SEXP Svalue);
SEXP RlpSolve_get_mat(SEXP Slp, SEXP Srow, SEXP Scolumn);
SEXP RlpSolve_set_obj_bound(SEXP Slp, SEXP Sobj_bound);
SEXP RlpSolve_get_obj_bound(SEXP Slp);
/*set_obj_fn*/
SEXP RlpSolve_set_obj_fnex(SEXP Slp, SEXP Srow, SEXP Scolno);
/*str_set_obj_fn*/
/*set_obj*/
SEXP RlpSolve_set_rh(SEXP Slp, SEXP Srow, SEXP Svalue);
SEXP RlpSolve_get_rh(SEXP Slp, SEXP Srow);
SEXP RlpSolve_set_rh_range(SEXP Slp, SEXP Srows, SEXP Sdeltavalue);
SEXP RlpSolve_get_rh_range(SEXP Slp, SEXP Srows);
SEXP RlpSolve_set_rh_vec(SEXP Slp, SEXP Srh);
/*str_set_rh_vec*/
SEXP RlpSolve_set_row_names(SEXP Slp, SEXP Srows, SEXP Snew_names);
SEXP RlpSolve_get_row_names(SEXP Slp, SEXP Srows);
SEXP RlpSolve_get_origrow_names(SEXP Slp, SEXP Srows);
SEXP RlpSolve_set_semicont(SEXP Slp, SEXP Scolumns, SEXP Ssc);
SEXP RlpSolve_is_semicont(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_upbo(SEXP Slp, SEXP Scolumns, SEXP Svalues);
SEXP RlpSolve_get_upbo(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_var_branch(SEXP Slp, SEXP Scolumns, SEXP Sbranch_mode);
SEXP RlpSolve_get_var_branch(SEXP Slp, SEXP Scolumns);
SEXP RlpSolve_set_var_weights(SEXP Slp, SEXP Sweights);


/*******************************
  * Solver settings
*******************************/

SEXP RlpSolve_default_basis(SEXP Slp);
/*read_basis*/
SEXP RlpSolve_reset_basis(SEXP Slp);
/*write_basis*/
SEXP RlpSolve_guess_basis(SEXP SLP, SEXP Sguessvector);
/*read_params*/
/*write_params*/
SEXP RlpSolve_reset_params(SEXP Slp);
SEXP RlpSolve_set_anti_degen(SEXP Slp, SEXP Santi_degen);
SEXP RlpSolve_is_anti_degen(SEXP Slp, SEXP Stestmasks);
SEXP RlpSolve_set_basis(SEXP Slp, SEXP Sbascolumn, SEXP Snonbasic);
SEXP RlpSolve_get_basis(SEXP Slp, SEXP Snonbasic);
SEXP RlpSolve_set_basiscrash(SEXP Slp, SEXP Smode);
SEXP RlpSolve_get_basiscrash(SEXP Slp);
SEXP RlpSolve_set_bb_depthlimit(SEXP Slp, SEXP Sbb_maxlevel);
SEXP RlpSolve_get_bb_depthlimit(SEXP Slp);
SEXP RlpSolve_set_bb_floorfirst(SEXP Slp, SEXP Sbb_floorfirst);
SEXP RlpSolve_get_bb_floorfirst(SEXP Slp);
SEXP RlpSolve_set_bb_rule(SEXP Slp, SEXP Sbb_rule);
SEXP RlpSolve_get_bb_rule(SEXP Slp);
/*set_BFP*/
/*has_BFP*/
/*is_nativeBFP*/
SEXP RlpSolve_set_break_at_first(SEXP Slp, SEXP Sbreak_at_first);
SEXP RlpSolve_is_break_at_first(SEXP Slp);
SEXP RlpSolve_set_break_at_value(SEXP Slp, SEXP Sbreak_at_value);
SEXP RlpSolve_get_break_at_value(SEXP Slp);
SEXP RlpSolve_set_epsb(SEXP Slp, SEXP Sepsb);
SEXP RlpSolve_get_epsb(SEXP Slp);
SEXP RlpSolve_set_epsd(SEXP Slp, SEXP Sepsd);
SEXP RlpSolve_get_epsd(SEXP Slp);
SEXP RlpSolve_set_epsel(SEXP Slp, SEXP Sepsel);
SEXP RlpSolve_get_epsel(SEXP Slp);
SEXP RlpSolve_set_epsint(SEXP Slp, SEXP Sepsint);
SEXP RlpSolve_get_epsint(SEXP Slp);
SEXP RlpSolve_set_epsperturb(SEXP Slp, SEXP Sepsperturb);
SEXP RlpSolve_get_epsperturb(SEXP Slp);
SEXP RlpSolve_set_epspivot(SEXP Slp, SEXP Sepspivot);
SEXP RlpSolve_get_epspivot(SEXP Slp);
SEXP RlpSolve_set_epslevel(SEXP Slp, SEXP Sepslevel);
SEXP RlpSolve_set_improve(SEXP Slp, SEXP Simprove);
SEXP RlpSolve_get_improve(SEXP Slp);
SEXP RlpSolve_set_maxim(SEXP Slp);
SEXP RlpSolve_is_maxim(SEXP Slp);
SEXP RlpSolve_set_maxpivot(SEXP Slp, SEXP Smax_num_inv);
SEXP RlpSolve_get_maxpivot(SEXP Slp);
SEXP RlpSolve_set_minim(SEXP Slp);
SEXP RlpSolve_set_mip_gap(SEXP Slp, SEXP Sabsolute, SEXP Smip_gap);
SEXP RlpSolve_get_mip_gap(SEXP Slp, SEXP Sabsolute);
SEXP RlpSolve_set_negrange(SEXP Slp, SEXP Snegrange);
SEXP RlpSolve_get_negrange(SEXP Slp);
SEXP RlpSolve_set_obj_in_basis(SEXP Slp, SEXP Sobj_in_basis);
SEXP RlpSolve_is_obj_in_basis(SEXP Slp);
SEXP RlpSolve_set_pivoting(SEXP Slp, SEXP Spivoting);
SEXP RlpSolve_get_pivoting(SEXP Slp);
SEXP RlpSolve_is_piv_mode(SEXP Slp, SEXP Stestmasks);
SEXP RlpSolve_is_piv_rule(SEXP Slp, SEXP Srules);
SEXP RlpSolve_set_preferdual(SEXP Slp, SEXP Sdodual);
SEXP RlpSolve_set_presolve(SEXP Slp, SEXP Sdo_presolve, SEXP Smaxloops);
SEXP RlpSolve_get_presolve(SEXP Slp);
SEXP RlpSolve_get_presolveloops(SEXP Slp);
SEXP RlpSolve_is_presolve(SEXP Slp, SEXP Stestmasks);
SEXP RlpSolve_set_scalelimit(SEXP Slp, SEXP Sscalelimit);
SEXP RlpSolve_get_scalelimit(SEXP Slp);
SEXP RlpSolve_set_scaling(SEXP Slp, SEXP Sscalemode);
SEXP RlpSolve_get_scaling(SEXP Slp);
SEXP RlpSolve_is_integerscaling(SEXP Slp);
SEXP RlpSolve_is_scalemode(SEXP Slp, SEXP Stestmasks);
SEXP RlpSolve_is_scaletype(SEXP Slp, SEXP Sscaletype);
SEXP RlpSolve_set_sense(SEXP Slp, SEXP Smaximize);
SEXP RlpSolve_set_simplextype(SEXP Slp, SEXP Ssimplextype);
SEXP RlpSolve_get_simplextype(SEXP Slp);
SEXP RlpSolve_set_solutionlimit(SEXP Slp, SEXP Slimit);
SEXP RlpSolve_get_solutionlimit(SEXP Slp);
SEXP RlpSolve_set_timeout(SEXP Slp, SEXP Ssectimeout);
SEXP RlpSolve_get_timeout(SEXP Slp);
SEXP RlpSolve_set_use_names(SEXP Slp, SEXP Sisrow, SEXP Suse_names);
SEXP RlpSolve_is_use_names(SEXP Slp, SEXP Sisrow);
SEXP RlpSolve_unscale(SEXP Slp);


/*******************************
  * Solve
*******************************/

SEXP RlpSolve_solve(SEXP Slp);
/*SEXP RlpSolve_lag_solve(SEXP Slp, SEXP Sstart_bound, SEXP Snum_iter);*/


/*******************************
  * Solution
*******************************/

SEXP RlpSolve_get_constraints(SEXP Slp);
/*get_ptr_constraints*/
/*get_constr_value*/
SEXP RlpSolve_get_objective(SEXP Slp);
SEXP RlpSolve_get_primal_solution(SEXP Slp);
/*get_ptr_primal_solution*/
SEXP RlpSolve_get_var_primalresult(SEXP Slp);
SEXP RlpSolve_get_sensitivity_obj(SEXP Slp);
/*get_ptr_sensitivity_obj*/
SEXP RlpSolve_get_sensitivity_objex(SEXP Slp);
/*get_ptr_sensitivity_objex*/
SEXP RlpSolve_get_sensitivity_rhs(SEXP Slp);
/*get_ptr_sensitivity_rhs*/
SEXP RlpSolve_get_dual_solution(SEXP Slp);
/*get_ptr_dual_solution*/
/*get_var_dualresult*/
SEXP RlpSolve_get_solutioncount(SEXP Slp);
SEXP RlpSolve_get_total_iter(SEXP Slp);
SEXP RlpSolve_get_total_nodes(SEXP Slp);
SEXP RlpSolve_get_variables(SEXP Slp);
/*get_ptr_variables*/
/*get_working_objective*/
/*is_feasible*/


/*******************************
  * Debug/print settings
*******************************/

/*set_debug*/
/*is_debug*/
/*set_lag_trace*/
/*is_lag_trace*/
/*set_outputstream*/
/*set_outputfile*/
/*set_print_sol*/
/*get_print_sol*/
/*set_trace*/
/*is_trace*/
SEXP RlpSolve_set_verbose(SEXP Slp, SEXP Sverbose);
SEXP RlpSolve_get_verbose(SEXP Slp);


/*******************************
  * Debug/print
*******************************/

/*print_constraints*/
/*print_debugdump*/
/*print_duals*/
/*print_lp*/
/*print_objective*/
/*print_scales*/
/*print_solution*/
/*print_str*/
/*SEXP RlpSolve_print_tableau(SEXP Slp);*/


/*******************************
  * Write model to file
*******************************/

SEXP RlpSolve_write_lp(SEXP Slp, SEXP Sfilename);
/*write_LP*/
/*write_lpex*/
SEXP RlpSolve_write_mps(SEXP Slp, SEXP Sfilename);
SEXP RlpSolve_write_freemps(SEXP Slp, SEXP Sfilename);
/*write_MPS*/
/*write_freeMPS*/
/*MPS_writefileex*/
/*write_XLI*/
/*set_XLI*/
/*has_XLI*/
/*is_nativeXLI*/

/*******************************
  * Miscellaneous routines
*******************************/

/*column_in_lp*/
SEXP RlpSolve_dualize_lp(SEXP Slp);
SEXP RlpSolve_get_lp_index(SEXP Slp, SEXP Sorig_indices);
/*SEXP RlpSolve_get_Lrows(SEXP Slp);*/
SEXP RlpSolve_get_Ncolumns(SEXP Slp);
SEXP RlpSolve_get_nonzeros(SEXP Slp);
SEXP RlpSolve_get_Norig_columns(SEXP Slp);
SEXP RlpSolve_get_Norig_rows(SEXP Slp);
SEXP RlpSolve_get_Nrows(SEXP Slp);
SEXP RlpSolve_get_orig_index(SEXP Slp, SEXP Slp_indices);
SEXP RlpSolve_get_status(SEXP Slp);
SEXP RlpSolve_get_statustext(SEXP Slp, SEXP Sstatuscode);
SEXP RlpSolve_lp_solve_version();
/*set_basisvar*/
SEXP RlpSolve_time_elapsed(SEXP Slp);




