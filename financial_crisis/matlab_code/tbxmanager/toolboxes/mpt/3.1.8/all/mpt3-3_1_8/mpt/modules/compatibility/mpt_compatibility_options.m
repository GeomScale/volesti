function S = mpt_compatibility_options
%
% Option settings for the compatibility module.
%

S.qpsolver = 1; % quadprog
S.lpsolver = 3; % cdd
S.milpsolver = -1;

% tolerances
S.abs_tol = 1e-7;
S.rel_tol = 1e-6;
S.step_size = 1e-5;
S.infbox = 10000;

% set to 'true' if you experience problems with solve QPs
S.use_old_solveQP = false;

% set to 'true' if you experience problems with solve LPs
S.use_old_solveLP = false;

% debugging/verbosity
S.debug_level = 1;
S.verbose = 1;
S.rescueQP = 0;

end
